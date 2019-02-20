#ifndef _PMS_H_
#define _PMS_H_

#include "basis_pms.h"
#include "deci.h"
#include <sstream>
Satlike::Satlike()
{
	problem_weighted=1;
	partial=1; //1 if the instance has hard clauses, and 0 otherwise.

	max_clause_length=0;
	min_clause_length=100000000;

	//size of the instance
	num_vars=0;		//var index from 1 to num_vars
	num_clauses=0;		//clause index from 0 to num_clauses-1
	num_hclauses=0;
	num_sclauses=0;

	print_time=240;
	cutoff_time=301;
	
	print1=0;
	print2=0;
	print_time1=50;
	print_time2=260;
	cutoff_time1=60;
	cutoff_time2=300;
}

void Satlike::settings()
{
	//steps
	weight_multiplier=100;
	tries=0;
	max_tries = 100000000;
	max_flips = 200000000;
	max_non_improve_flip = 10000000;
	max_flips = max_non_improve_flip;
	step=1;

	large_clause_count_threshold=0;
	soft_large_clause_count_threshold=0;

	rdprob=0.01;
	hd_count_threshold=42;
	rwprob=0.091; 
	smooth_probability=0.000003;
	h_inc=1;
	softclause_weight_threshold=400;

    cfg_bayesian_init_epochs = 10;
	
	if(top_clause_weight==555 || top_clause_weight==3345){max_non_improve_flip = 100000000;max_flips = max_non_improve_flip;return;}//for tpr-multiple-path
	if(num_vars==21368 && max_clause_length==8)//for aes
	{
		//weight_multiplier=1;
		max_non_improve_flip = 100000000;
		max_flips = max_non_improve_flip;
		h_inc=2;
		softclause_weight_threshold=0;
		rdprob=0.01;
		hd_count_threshold=15;
		rwprob=0; 
		smooth_probability=0.0001;
		return ;
	}
	if(num_vars>12000 && num_vars<17000 && top_clause_weight==49771 && max_clause_length>310 && max_clause_length<351)//for fault
	{
		//weight_multiplier=1;
		h_inc=2;
		softclause_weight_threshold=0;
		rdprob=0.01;
		hd_count_threshold=15;
		rwprob=0; 
		smooth_probability=0.0001;
		return ;
	}

	if(num_vars<1100) //level 2
	{
		weight_multiplier=1;
		h_inc=2;
		softclause_weight_threshold=0;
		rdprob=0.01;
		hd_count_threshold=15;
		rwprob=0; 
		smooth_probability=0.01;
		return ;
	}
	
}

void Satlike::allocate_memory()
{
	int malloc_var_length = num_vars+10;
	int malloc_clause_length = num_clauses+10;
	
	unit_clause = new lit[malloc_clause_length];

	var_lit = new lit* [malloc_var_length];
	var_lit_count = new int [malloc_var_length];
	clause_lit = new lit* [malloc_clause_length];
	clause_lit_count = new int [malloc_clause_length];
	
	score = new long long [malloc_var_length];
	var_neighbor = new int* [malloc_var_length];
	var_neighbor_count = new int [malloc_var_length];
	time_stamp = new long long [malloc_var_length];
	neighbor_flag = new int [malloc_var_length];
	temp_neighbor = new int [malloc_var_length];
	
	org_clause_weight = new long long [malloc_clause_length];
	clause_weight = new long long [malloc_clause_length];
	sat_count = new int [malloc_clause_length];
	sat_var = new int [malloc_clause_length];
	clause_selected_count = new long long [malloc_clause_length];
	best_soft_clause = new int [malloc_clause_length];
	
	hardunsat_stack = new int [malloc_clause_length];
	index_in_hardunsat_stack = new int [malloc_clause_length];
	softunsat_stack = new int [malloc_clause_length];
	index_in_softunsat_stack = new int [malloc_clause_length];
	
	unsatvar_stack = new int [malloc_var_length];
	index_in_unsatvar_stack = new int [malloc_var_length];
	unsat_app_count = new int [malloc_var_length];
	
	goodvar_stack = new int [malloc_var_length];
	already_in_goodvar_stack = new int[malloc_var_length];

	cur_soln = new int [malloc_var_length];
	best_soln = new int [malloc_var_length];
	local_opt_soln = new int[malloc_var_length];
	
	large_weight_clauses = new int [malloc_clause_length];
	soft_large_weight_clauses = new int [malloc_clause_length];
	already_in_soft_large_weight_stack = new int [malloc_clause_length];

	best_array = new int [malloc_var_length];
	temp_lit = new int [malloc_var_length];

    parameters = new BetaDist[malloc_var_length];
    updated_params = new BetaDist[malloc_var_length];
}

void Satlike::free_memory()
{
	int i;
	for (i = 0; i < num_clauses; i++) 
		delete[] clause_lit[i];
	
	for(i=1; i<=num_vars; ++i)
	{
		delete[] var_lit[i];
		delete[] var_neighbor[i];
	}
	
	delete [] var_lit;
	delete [] var_lit_count;
	delete [] clause_lit;
	delete [] clause_lit_count;

	delete [] score;
	delete [] var_neighbor;
	delete [] var_neighbor_count;
	delete [] time_stamp;
	delete [] neighbor_flag;
	delete [] temp_neighbor;
	
	delete [] org_clause_weight;
	delete [] clause_weight;
	delete [] sat_count;
	delete [] sat_var;
	delete [] clause_selected_count;
	delete [] best_soft_clause;
	
	delete [] hardunsat_stack;
	delete [] index_in_hardunsat_stack;
	delete [] softunsat_stack;
	delete [] index_in_softunsat_stack;
	
	delete [] unsatvar_stack;
	delete [] index_in_unsatvar_stack;
	delete [] unsat_app_count;
	
	delete [] goodvar_stack;
	delete [] already_in_goodvar_stack;
	
	//delete [] fix;
	delete [] cur_soln;
	delete [] best_soln;
	delete [] local_opt_soln;
	
	delete [] large_weight_clauses;
	delete [] soft_large_weight_clauses;
	delete [] already_in_soft_large_weight_stack;

	delete [] best_array;
	delete [] temp_lit;

    delete [] parameters;
    delete [] updated_params;
}

void Satlike::build_neighbor_relation()
{
	int	i,j,count;
	int v,c,n;
    int temp_neighbor_count;

	for(v=1; v<=num_vars; ++v)
	{
		neighbor_flag[v] = 1;
		temp_neighbor_count = 0;
		
		for(i=0; i<var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
			for(j=0; j<clause_lit_count[c]; ++j)
			{
                n=clause_lit[c][j].var_num;
				if(neighbor_flag[n]!=1)
				{
					neighbor_flag[n] = 1;
					temp_neighbor[temp_neighbor_count++] = n;
				}
			}
		}
        
		neighbor_flag[v] = 0;
        
		var_neighbor[v] = new int[temp_neighbor_count];
		var_neighbor_count[v]=temp_neighbor_count;
        
		count = 0;
		for(i=0; i<temp_neighbor_count; i++)
		{
			var_neighbor[v][count++] = temp_neighbor[i];
			neighbor_flag[temp_neighbor[i]] = 0;
		}
	}
}

void Satlike::build_instance(char *filename)
{
	char    line[1024];
	string  line2;
	char    tempstr1[10];
	char    tempstr2[10];
	int     cur_lit;
	int     i,v,c;
    //int     temp_lit[MAX_VARS];
	
	ifstream infile(filename);
	if(infile==NULL) 
    {
		cout<<"c the input filename "<<filename<<" is invalid, please input the correct filename."<<endl;
		exit(-1);
	}
    
	/*** build problem data structures of the instance ***/
	
	getline(infile,line2);
	
	while(line2[0]!='p')
	{
		getline(infile,line2);
	}
	for(i=0;i<1024;i++)
	{
		line[i]=line2[i];
	}
    int read_items;
    
	read_items=sscanf(line, "%s %s %d %d %lld", tempstr1, tempstr2, &num_vars, &num_clauses, &top_clause_weight);
    
    if(read_items<5)
    {
    	if(strcmp(tempstr2,"cnf")==0) 
		{			
			problem_weighted=0;
		}
	
		partial = 0;
		top_clause_weight=-1;
	}

	allocate_memory();

	for (c = 0; c < num_clauses; c++) 
	{
		clause_lit_count[c] = 0;
		clause_lit[i]=NULL;
	}
	for (v=1; v<=num_vars; ++v)
	{
		var_lit_count[v] = 0;
		var_lit[v]=NULL;
		var_neighbor[v]=NULL;
	}
    
    num_hclauses=num_sclauses=0;
	unit_clause_count=0;
	//Now, read the clauses, one at a time.
	int lit_redundent,clause_redundent;
	for (c = 0; c < num_clauses; ) 
	{
		clause_redundent=0;
        clause_lit_count[c] = 0;
        
        if (problem_weighted!=0) 
        	infile>>org_clause_weight[c];
        else org_clause_weight[c]=1;
        
        if (org_clause_weight[c]!=top_clause_weight) 
        {
            total_soft_weight += org_clause_weight[c];
            num_sclauses++;
        }
        else num_hclauses++;
		
		infile>>cur_lit;
		while (cur_lit != 0) { 
            
			lit_redundent=0;
			for(int p=0; p<clause_lit_count[c]; p++)
			{
				if(cur_lit==temp_lit[p]){
					lit_redundent=1;
					break;
				}
				else if(cur_lit==-temp_lit[p]){
					clause_redundent=1;
					break;
				}
			}
			
			if(lit_redundent==0)
			{
				temp_lit[clause_lit_count[c]] = cur_lit;
				clause_lit_count[c]++;
			}
            
			infile>>cur_lit;
		}
		
		if(clause_redundent==0) //the clause is not tautology
		{
			clause_lit[c] = new lit[clause_lit_count[c]+1];
            
			for(i=0; i<clause_lit_count[c]; ++i)
			{
				clause_lit[c][i].clause_num = c;
				clause_lit[c][i].var_num = abs(temp_lit[i]);
				
				if (temp_lit[i] > 0) clause_lit[c][i].sense = 1;
                else clause_lit[c][i].sense = 0;
                
				var_lit_count[clause_lit[c][i].var_num]++;
			}
            clause_lit[c][i].var_num=0; 
            clause_lit[c][i].clause_num = -1;
            
			if(clause_lit_count[c]==1)
			{
				unit_clause[unit_clause_count++]=clause_lit[c][0];
			}

            if(clause_lit_count[c]>max_clause_length) max_clause_length=clause_lit_count[c];
			if(clause_lit_count[c]<min_clause_length) min_clause_length=clause_lit_count[c];
            
			c++;
		}
		else 
		{
			num_clauses--;
			clause_lit_count[c] = 0;
		}
	}
	infile.close();
	//creat var literal arrays
	for (v=1; v<=num_vars; ++v)
	{
		var_lit[v] = new lit[var_lit_count[v]+1];
		var_lit_count[v] = 0;	//reset to 0, for build up the array
	}
	//scan all clauses to build up var literal arrays
	for (c = 0; c < num_clauses; ++c) 
	{
		for(i=0; i<clause_lit_count[c]; ++i)
		{
			v = clause_lit[c][i].var_num;
			var_lit[v][var_lit_count[v]] = clause_lit[c][i];
			++var_lit_count[v];
		}
	}
    for (v=1; v<=num_vars; ++v)
		var_lit[v][var_lit_count[v]].clause_num=-1;
        

	build_neighbor_relation();

    best_soln_feasible=0;
    opt_unsat_weight=total_soft_weight+1;

    bayesian();

    for( int v=1; v<=num_vars; v++ )
        best_soln[v] = (parameters[v].a > parameters[v].b);
}

void Satlike::bayesian_update(int c)
{
    double coeff_product = 1.0;
    for( int j=0; j<clause_lit_count[c]; j++ )
    {
        int v = clause_lit[c][j].var_num;
        int sgn = clause_lit[c][j].sense;

        updated_params[v].a = parameters[v].a + (!sgn);
        updated_params[v].b = parameters[v].b + (sgn);

        double coeff = (sgn ? parameters[v].b : parameters[v].a) / (parameters[v].a + parameters[v].b);
        coeff_product *= coeff;
    }
    double normalization_constant = 1 - coeff_product;

    for( int j=0; j<clause_lit_count[c]; j++ )
    {
        int v = clause_lit[c][j].var_num;

        double sumP = parameters[v].a + parameters[v].b;
        double sumUP = updated_params[v].a + updated_params[v].b;

        double *p[2], *up[2];
        p[0] = &parameters[v].a;
        p[1] = &parameters[v].b;
        up[0] = &updated_params[v].a;
        up[1] = &updated_params[v].b;
        for( int k=0; k<=1; k++ )
        {
            double moment1 = *p[k] / sumP - coeff_product * *up[k] / sumUP;
            double moment2 = *p[k] * (*p[k] + 1) / ((sumP) * (sumP+1)) - coeff_product * *up[k] * (*up[k] + 1) / ((sumUP) * (sumUP+1));

            moment1 /= normalization_constant;
            moment2 /= normalization_constant;

            *p[k] = moment1 * (moment1 - moment2) / (moment2 - moment1*moment1);
        }
    }
}

void Satlike::bayesian()
{
    int epochs = cfg_bayesian_init_epochs;
    int n = num_vars;

    for( int i=0; i<n; i++ )
    {
        parameters[i].a = 1 + (double)rand() / RAND_MAX;
        parameters[i].b = 1 + (double)rand() / RAND_MAX;
    }

    for( int k=0; k<epochs; k++ )
    {
        for( int i=0; i<num_clauses; i++ )
            bayesian_update(i);
    }
}


void Satlike::init(vector<int>& init_solution)
{
	int 		v,c;
	int			i,j;
    hard_unsat_nb=0;
    soft_unsat_weight=0;
    local_soln_feasible=0;
	local_opt_unsat_weight=top_clause_weight+1;

	if(best_soln_feasible==0)
	{
		for (c = 0; c<num_clauses; c++)
    	{
        already_in_soft_large_weight_stack[c]=0;
		clause_selected_count[c]=0;

		if (org_clause_weight[c]==top_clause_weight) 
            clause_weight[c] = 1;
        else
		{
			clause_weight[c] = 0;
		}
    	}
	}
	else
	{
    //Initialize clause information
	for (c = 0; c<num_clauses; c++)
    {
        already_in_soft_large_weight_stack[c]=0;
		clause_selected_count[c]=0;

		if (org_clause_weight[c]==top_clause_weight) 
            clause_weight[c] = 1;
        else
		{
			clause_weight[c] = org_clause_weight[c]*weight_multiplier;;
			if(clause_weight[c]>1 && already_in_soft_large_weight_stack[c]==0)
			{
				already_in_soft_large_weight_stack[c]=1;
				soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
			}
		}
    }
    }

    //init solution
    if(init_solution.size()==0)
    {
        for (v = 1; v <= num_vars; v++) 
        {
            cur_soln[v]=0;
		    time_stamp[v] = 0;
		    unsat_app_count[v] = 0;
	    }
    }
	else 
    {
        for (v = 1; v <= num_vars; v++) 
        {
            cur_soln[v]=init_solution[v];
		    time_stamp[v] = 0;
		    unsat_app_count[v] = 0;
	    }
    }

	//init stacks
	hardunsat_stack_fill_pointer = 0;
    softunsat_stack_fill_pointer = 0;
	unsatvar_stack_fill_pointer = 0;
	large_weight_clauses_count=0;
    soft_large_weight_clauses_count=0;

	/* figure out sat_count, sat_var and init unsat_stack */
	for (c=0; c<num_clauses; ++c) 
	{
		sat_count[c] = 0;
		for(j=0; j<clause_lit_count[c]; ++j)
		{
			if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
			{
				sat_count[c]++;
				sat_var[c] = clause_lit[c][j].var_num;	
			}
		}
		if (sat_count[c] == 0) 
        {
            unsat(c);
        }
	}

	/*figure out score*/
	for (v=1; v<=num_vars; v++) 
	{
		score[v]=0;
		for(i=0; i<var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
            if (sat_count[c]==0) score[v]+=clause_weight[c];
            else if (sat_count[c]==1 && var_lit[v][i].sense==cur_soln[v]) score[v]-=clause_weight[c];
		}
	}

	//init goodvars stack
	goodvar_stack_fill_pointer = 0;
	for (v=1; v<=num_vars; v++) 
	{
		if(score[v]>0)
		{
			already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
			mypush(v,goodvar_stack);
		}
		else
			already_in_goodvar_stack[v] = -1;
	}
}

int Satlike::pick_var()
{
	int     i,v;
	int     best_var;
    
	if(goodvar_stack_fill_pointer>0 )
	{
		if( (rand()%MY_RAND_MAX_INT)*BASIC_SCALE< rdprob ) 
            return goodvar_stack[rand()%goodvar_stack_fill_pointer]; 
            
        if (goodvar_stack_fill_pointer < hd_count_threshold)
        {
		    best_var = goodvar_stack[0];
		    for (i=1; i<goodvar_stack_fill_pointer; ++i) 
		    {
		    	v = goodvar_stack[i];
		    	if (score[v]>score[best_var]) best_var=v;
		        else if (score[v]==score[best_var])
		        {   
		            if(time_stamp[v]<time_stamp[best_var]) best_var=v;
		        } 
		    }
		    return best_var;
        }
        else
		{
		    best_var = goodvar_stack[rand()%goodvar_stack_fill_pointer];
			for (i=1; i<hd_count_threshold; ++i) 
		    {
		        v = goodvar_stack[rand()%goodvar_stack_fill_pointer];
		        if (score[v]>score[best_var]) best_var=v;
		        else if (score[v]==score[best_var])
		        {   
		        	if(time_stamp[v]<time_stamp[best_var]) best_var=v;
		        } 
		    }
		    return best_var;
        }
	}
	
    update_clause_weights();
    
    int sel_c;    
    lit* p;
    
    if (hardunsat_stack_fill_pointer>0) 
    {
    	sel_c = hardunsat_stack[rand()%hardunsat_stack_fill_pointer]; 
    }
    else 
    {
    	sel_c= softunsat_stack[rand()%softunsat_stack_fill_pointer];

		long long best_selected_num;
		int best_soft_clause_count;
    	int c;
    	best_soft_clause[0]=softunsat_stack[0];
    	best_selected_num=clause_selected_count[softunsat_stack[0]];
    	best_soft_clause_count=1;
    	
    	for(int j=1;j<softunsat_stack_fill_pointer;++j)
    	{
    		c=softunsat_stack[j];
    		if(clause_selected_count[c]<best_selected_num)
    		{
    			best_soft_clause[0]=c;
    			best_selected_num=clause_selected_count[c];
    			best_soft_clause_count=1;
    		}
    		else if(clause_selected_count[c]==best_selected_num)
    		{
    			best_soft_clause[best_soft_clause_count]=c;
    			best_soft_clause_count++;
    		}
    	}
    	sel_c=best_soft_clause[rand()%best_soft_clause_count];
    	clause_selected_count[sel_c]++;
    }
    if( (rand()%MY_RAND_MAX_INT)*BASIC_SCALE< rwprob )  
        return clause_lit[sel_c][rand()%clause_lit_count[sel_c]].var_num;
 
    best_var = clause_lit[sel_c][0].var_num;
    p=clause_lit[sel_c];
	for(p++; (v=p->var_num)!=0; p++)  
    {           
        if (score[v]>score[best_var]) best_var=v;
		else if (score[v]==score[best_var])
		{   
			if(time_stamp[v]<time_stamp[best_var]) best_var=v;
		} 
    }
    return best_var;
}

void Satlike::update_goodvarstack1(int flipvar)
{
	int v;
	//remove the vars no longer goodvar in goodvar stack 
	for(int index=goodvar_stack_fill_pointer-1; index>=0; index--)
	{
		v = goodvar_stack[index];
		if(score[v]<=0)
		{
			goodvar_stack[index] = mypop(goodvar_stack);
			already_in_goodvar_stack[v] = -1;
		}	
	}

	//add goodvar
	for(int i=0; i<var_neighbor_count[flipvar]; ++i)
	{
		v = var_neighbor[flipvar][i];
		if( score[v] > 0) 
		{
            if(already_in_goodvar_stack[v] == -1)
            {
            	already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v,goodvar_stack);
            }
		}
	}
}
void Satlike::update_goodvarstack2(int flipvar)
{
	if(score[flipvar]>0 && already_in_goodvar_stack[flipvar]==-1)
	{
		already_in_goodvar_stack[flipvar]=goodvar_stack_fill_pointer;
		mypush(flipvar,goodvar_stack);
	}
	else if(score[flipvar]<=0 && already_in_goodvar_stack[flipvar]!=-1)
	{
		int index=already_in_goodvar_stack[flipvar];
		int last_v=mypop(goodvar_stack);
		goodvar_stack[index]=last_v;
		already_in_goodvar_stack[last_v]=index;
		already_in_goodvar_stack[flipvar]=-1;	
	}
	int i,v;
	for(i=0; i<var_neighbor_count[flipvar]; ++i)
	{
		v = var_neighbor[flipvar][i];
		if( score[v] > 0) 
		{
            if(already_in_goodvar_stack[v] == -1)
            {
            	already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v,goodvar_stack);
            }
		}
		else if(already_in_goodvar_stack[v]!=-1)
		{
			int index=already_in_goodvar_stack[v];
			int last_v=mypop(goodvar_stack);
			goodvar_stack[index]=last_v;
			already_in_goodvar_stack[last_v]=index;
			already_in_goodvar_stack[v]=-1;	
		}
	}
}

void Satlike::flip(int flipvar)
{
	int i,v,c;
	int index;
	lit* clause_c;

    int	org_flipvar_score = score[flipvar];
	cur_soln[flipvar] = 1 - cur_soln[flipvar];
	
    for(i=0; i<var_lit_count[flipvar]; ++i)
	{
		c = var_lit[flipvar][i].clause_num;
		clause_c = clause_lit[c];
        
        if(cur_soln[flipvar] == var_lit[flipvar][i].sense)
        {
            ++sat_count[c];
            if (sat_count[c] == 2) //sat_count from 1 to 2
            {
                score[sat_var[c]] += clause_weight[c];
            }
            else if (sat_count[c] == 1) // sat_count from 0 to 1
            {
                sat_var[c] = flipvar;//record the only true lit's var
                for(lit* p=clause_c; (v=p->var_num)!=0; p++) 
                {
                	score[v] -= clause_weight[c];
                }
                sat(c);
            }
        }
        else // cur_soln[flipvar] != cur_lit.sense
        {
            --sat_count[c];
            if (sat_count[c] == 1) //sat_count from 2 to 1
            {
                for(lit* p=clause_c; (v=p->var_num)!=0; p++) 
                {
                     if(p->sense == cur_soln[v] )
                    {
                        score[v] -= clause_weight[c];
                        sat_var[c] = v;
                        break;
                    }
                }
            }
            else if (sat_count[c] == 0) //sat_count from 1 to 0
            {
                for(lit* p=clause_c; (v=p->var_num)!=0; p++)
                {
					score[v] += clause_weight[c];
				}    
                unsat(c);    
            }//end else if  
        }//end else
	}

	//update information of flipvar
	score[flipvar] = -org_flipvar_score;
	update_goodvarstack1(flipvar);
}

void Satlike::local_search(vector<int>& init_solution)
{
	settings();
	max_flips = 200000000;
	for(tries=1;tries<max_tries;++tries)
	{
		init(init_solution);
		for (step = 1; step<max_flips; ++step)
		{
			if (hard_unsat_nb==0 && (soft_unsat_weight<opt_unsat_weight || best_soln_feasible==0) )
        	{
		    	if (soft_unsat_weight<top_clause_weight) 
		    	{
		        	best_soln_feasible = 1;
		        	opt_unsat_weight = soft_unsat_weight;
					opt_time = get_runtime(); 
		        	for(int v=1; v<=num_vars; ++v) best_soln[v] = cur_soln[v];
		    	}
        	}
			if (step%1000==0)
        	{
        		double elapse_time=get_runtime();
            	if(elapse_time>=cutoff_time) return;
				else if(opt_unsat_weight==0) return;    
        	}
			int flipvar = pick_var();
			flip(flipvar);
			time_stamp[flipvar] = step;
		}
	} 
}

void Satlike::print_best_solution()
{
	if (best_soln_feasible==0) return;

     /*
     printf("v");
     for (int i=1; i<=num_vars; i++) {
		printf(" ");
		if(best_soln[i]==0) printf("-");
		printf("%d", i);
     }
     printf("\n");
     */

     std::stringstream s;
     s << "v";

     for (int i=1; i<=num_vars; i++) {
       s << " ";
       if(best_soln[i]==0) s << "-";
       s << i;
     }
     printf("%s\n",s.str().c_str());
}

void Satlike::local_search_with_decimation(vector<int>& init_solution, char* inputfile)
{
	settings();
	
	Decimation deci(var_lit,var_lit_count,clause_lit,org_clause_weight,top_clause_weight);
	deci.make_space(num_clauses,num_vars);

	for(int i=1;i<=num_vars;++i)
	{
		local_opt_soln[i]=rand()%2;
	}

	for(tries=1;tries<max_tries;++tries)
	{
		deci.init(local_opt_soln,best_soln,unit_clause,unit_clause_count,clause_lit_count);

		deci.unit_prosess();

		init_solution.resize(num_vars+1);
		for(int i=1;i<=num_vars;++i)
		{
			init_solution[i]=deci.fix[i];
		}

		init(init_solution);
		for (step = 1; step<max_flips; ++step)
		{
			if (hard_unsat_nb==0 && (soft_unsat_weight<opt_unsat_weight || best_soln_feasible==0) )
        	{
        		if(best_soln_feasible==0)
				{
					best_soln_feasible=1;
					opt_unsat_weight = soft_unsat_weight;
					cout<<"o "<<opt_unsat_weight<<endl;
					//opt_time = get_runtime(); 
		        	for(int v=1; v<=num_vars; ++v) best_soln[v] = cur_soln[v];
		        	break;
				}
		    	if (soft_unsat_weight<top_clause_weight) 
		    	{
		        	
		        	opt_unsat_weight = soft_unsat_weight;
		        	cout<<"o "<<opt_unsat_weight<<endl;
		        	
					opt_time = get_runtime(); 
		        	for(int v=1; v<=num_vars; ++v) best_soln[v] = cur_soln[v];

				/*
		        	if(print1==1 || print2==1)
		        	{
		        		print_best_solution();
		        	}
				*/
		        	
		        	if(opt_unsat_weight==0)
		        	{
		        		print_best_solution();
		        		return;
		        	}
		    	}
        	}
        	
        	if(hard_unsat_nb==0 && (soft_unsat_weight<local_opt_unsat_weight || local_soln_feasible==0))
			{
				if(soft_unsat_weight<top_clause_weight)
				{
					local_soln_feasible=1;
					local_opt_unsat_weight=soft_unsat_weight;
					max_flips=step+max_non_improve_flip;
					for(int v=1;v<=num_vars;++v) local_opt_soln[v]=cur_soln[v];
				}
			}
			
			
			if (step%1000==0)
        	{
        		double elapse_time=get_runtime();
        		
        		if(print1==0 && elapse_time>print_time1 && elapse_time<60)
        		{
        			print1=1;
        			//if(best_soln_feasible==1)
        			//	print_best_solution();
        		}
        		
        		if(elapse_time > 60) print1=0;
        		
        		if(print2==0 && elapse_time>print_time2)
        		{
        			print2=1;
        			//if(best_soln_feasible==1)
        			//	print_best_solution();
        		}
        		
            	if(elapse_time>=cutoff_time) {deci.free_memory();return;}
        	}
        	
			int flipvar = pick_var();
			flip(flipvar);
			time_stamp[flipvar] = step;
		}
	} 
}

bool Satlike::verify_sol()
{        
	int c,j,flag;
    long long verify_unsat_weight=0;
    
	for (c = 0; c<num_clauses; ++c) 
	{
        flag = 0;
        for(j=0; j<clause_lit_count[c]; ++j)
            if (best_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense) {flag = 1; break;}
        
        if(flag ==0)
        {
            if(org_clause_weight[c]==top_clause_weight)//verify hard clauses
            {
                //output the clause unsatisfied by the solution
                cout<<"c Error: hard clause "<<c<<" is not satisfied"<<endl;

                cout<<"c ";
                for(j=0; j<clause_lit_count[c]; ++j)
                {
                    if(clause_lit[c][j].sense==0)cout<<"-";
                    cout<<clause_lit[c][j].var_num<<" ";
                }
                cout<<endl;
                cout<<"c ";
                for(j=0; j<clause_lit_count[c]; ++j)
                    cout<<best_soln[clause_lit[c][j].var_num]<<" ";
                cout<<endl;
                return 0;
            }
            else
            {
                verify_unsat_weight+=org_clause_weight[c];
            }
        }
    }

    if(verify_unsat_weight==opt_unsat_weight)  return 1;
    else{
    	cout<<"c Error: find opt="<<opt_unsat_weight<<", but verified opt="<<verify_unsat_weight<<endl;
    }
    return 0;
}

void Satlike::simple_print()
{
	if(best_soln_feasible==1)
	{
		if(verify_sol()==1)
			cout<<opt_unsat_weight<<'\t'<<opt_time<<endl;
		else
			cout<<"solution is wrong "<<endl;
	}
	else 
		cout<<-1<<'\t'<<-1<<endl;
}

void Satlike::increase_weights()
{
	int i,c,v;
    for(i=0; i < hardunsat_stack_fill_pointer; ++i)
    {
        c=hardunsat_stack[i];
		clause_weight[c]+=h_inc;

        if(clause_weight[c] == (h_inc+1))
			large_weight_clauses[large_weight_clauses_count++] = c;
		
        for(lit* p=clause_lit[c]; (v=p->var_num)!=0; p++)
		{
			score[v]+=h_inc;
            if (score[v]>0 && already_in_goodvar_stack[v]==-1)           
            {
                already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v,goodvar_stack);
            }
		}
    }
    if(best_soln_feasible!=0)
    {
	for(i=0; i < softunsat_stack_fill_pointer; ++i)
    {
        c=softunsat_stack[i];
        if(clause_weight[c]>softclause_weight_threshold) continue;
		else clause_weight[c]++;
		
        if(clause_weight[c]>1 && already_in_soft_large_weight_stack[c]==0)
        {
            already_in_soft_large_weight_stack[c]=1;
			soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
        }
        for(lit* p=clause_lit[c]; (v=p->var_num)!=0; p++)
		{
			score[v]++;
            if (score[v]>0 && already_in_goodvar_stack[v]==-1)           
            {
				already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v,goodvar_stack);
            }
		}
    }
    }
}

void Satlike::smooth_weights()
{
	int i, clause, v;

	for(i=0; i<large_weight_clauses_count; i++)
	{
		clause = large_weight_clauses[i];
		if(sat_count[clause]>0)
		{
			clause_weight[clause]-=h_inc;

			if(clause_weight[clause]==1)
			{
				large_weight_clauses[i] = large_weight_clauses[--large_weight_clauses_count];
				i--;
			}
			if(sat_count[clause] == 1)
			{
				v = sat_var[clause];
				score[v]+=h_inc;
                if (score[v]>0 && already_in_goodvar_stack[v]==-1)           
                {
					already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v,goodvar_stack);
                }
			}
		}
	}
	if(best_soln_feasible!=0)
	{
	for(i=0; i<soft_large_weight_clauses_count; i++)
	{
		clause = soft_large_weight_clauses[i];
		if(sat_count[clause]>0)
		{
			clause_weight[clause]--;
			if(clause_weight[clause]==1 && already_in_soft_large_weight_stack[clause]==1)
			{
                already_in_soft_large_weight_stack[clause]=0;
				soft_large_weight_clauses[i] = soft_large_weight_clauses[--soft_large_weight_clauses_count];
				i--;
			}
			if(sat_count[clause] == 1)
			{
				v = sat_var[clause];
				score[v]++;
                if (score[v]>0 && already_in_goodvar_stack[v]==-1)           
                {
					already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v,goodvar_stack);
                }
			}
		}
	}
	}
}

void Satlike::update_clause_weights()
{	
	if( ((rand()%MY_RAND_MAX_INT)*BASIC_SCALE)<smooth_probability 
       && large_weight_clauses_count>large_clause_count_threshold  )
	{
		smooth_weights();
	}
	else 
	{
		increase_weights();
	}
}

inline void Satlike::unsat(int clause)
{
    if(org_clause_weight[clause]==top_clause_weight) 
    {
        index_in_hardunsat_stack[clause] = hardunsat_stack_fill_pointer;
        mypush(clause,hardunsat_stack);
        hard_unsat_nb++;
    }
    else 
    {
        index_in_softunsat_stack[clause] = softunsat_stack_fill_pointer;
        mypush(clause,softunsat_stack);
        soft_unsat_weight += org_clause_weight[clause];
    }   
}

inline void Satlike::sat(int clause)
{
	int index,last_unsat_clause;

    if (org_clause_weight[clause]==top_clause_weight) 
    {
        last_unsat_clause = mypop(hardunsat_stack);
        index = index_in_hardunsat_stack[clause];
        hardunsat_stack[index] = last_unsat_clause;
        index_in_hardunsat_stack[last_unsat_clause] = index;

        hard_unsat_nb--;
    }
    else 
    {
        last_unsat_clause = mypop(softunsat_stack);
        index = index_in_softunsat_stack[clause];
        softunsat_stack[index] = last_unsat_clause;
        index_in_softunsat_stack[last_unsat_clause] = index;
        
        soft_unsat_weight -= org_clause_weight[clause];
    }
}

#endif
