#include "trans_news.h"
#include "mpi.h"
#include "string.h"

int myid,numprocs;
SubArea*myarea = NULL;
int**masterbuf1;int**masterbuf3;int**masterbuf5;
double**masterbuf2;double**masterbuf4;
int*infectious_id1;int*infectious_id2;
int*prevention_id1;int*prevention_id2;
double *infectious_time1;double* infectious_time2;
unordered_map<int,int>pid_slaveid; //pid-slave
ofstream ffout("result_trans.txt",ofstream::app);
void Initial_subarea();
void send_data_to_master(vector<pair<int,double> >& v,int);
void master_receive_from_slave();
void slave_receive_from_master();
int main(int argc,char**argv)
{
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	srand((myid+1)*(unsigned int)time(NULL));  //设置种子
	char* parafile = "iparameter.txt"; //file name
	ifstream parafin(parafile);
	Parameter(parafin);
	cout<<myid<<"ok"<<endl;
	int* seedsize;
	int** randseed;
	if(myid==0)
	{
		char* filename = ".//slave_pop.txt";
		ifstream fin(filename);
		int slaveid,pid,age;
		while(fin>>slaveid>>pid>>age)
		{
			pid_slaveid[pid] = slaveid;
		}
		fin.close();
		masterbuf1 = new int*[numprocs];
		masterbuf2 = new double*[numprocs];
		masterbuf3 = new int*[numprocs];
		masterbuf4 = new double*[numprocs];
		masterbuf5 = new int*[numprocs];
		for(int i=0;i<numprocs;i++)
		{
			masterbuf1[i] = new int[MAX_SLAVE];
			masterbuf2[i] = new double[MAX_SLAVE];
			masterbuf3[i] = new int[MAX_SLAVE];
			masterbuf4[i] = new double[MAX_SLAVE];
			masterbuf5[i] = new int[MAX_SLAVE];
		}
		cout<<"myid = "<<myid<<"slave_pop&masterbuf"<<endl;
	}
	else
	{
		myarea = new SubArea(myid);
		char* filename = ".//slave_pop.txt";
		ifstream fin(filename);
		int slaveid,pid,age;
		while(fin>>slaveid>>pid>>age)
		{
			if(slaveid == myid)
			{
				myarea->Individuals[pid] = new Individual(pid,age,slaveid);
			}
		}
		fin.close();
		readNet(".//JTNet.txt",1,myarea);
        readNet(".//XXNet.txt",2,myarea);
        readNet(".//GZNet.txt",3,myarea);
        readNet(".//PYNet.txt",4,myarea);
        readNet(".//SQNet.txt",5,myarea);
        readNet(".//TQNet.txt",6,myarea);
        infectious_id1 = new int[MAX_SLAVE];
        infectious_id2 = new int[MAX_SLAVE];
        prevention_id1 = new int[MAX_SLAVE];
        prevention_id2 = new int[MAX_SLAVE];
        infectious_time1 = new double[MAX_SLAVE];
        infectious_time2 = new double[MAX_SLAVE];
		cout<<"readNet"<<endl;
	} // first if myid!=0
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid ==0) //生成感染者
	{
		cout<<"hello"<<endl;
		vector<vector<int> > id;
		for(int i=0;i<numprocs;i++)
		{
			vector<int> which;
			id.push_back(which);
		}
		unordered_set<int> rnum;
        int pid_size = pid_slaveid.size();
		for(int i=0;i<SEED_NUM;i++)
		{
			int index = 0;
			while(1)
			{
				index = 1+(int)((double)rand()/((double)RAND_MAX+1.0)*(double)pid_size);
				if(rnum.find(index)!=rnum.end())
				{
					continue;
				}
				else
				{
					rnum.insert(index);
					break;
				}
			}
			id_age_iter it = pid_slaveid.find(index);
			if(it==pid_slaveid.end())
			{
				cout<<"can't find error"<<endl;
				exit(0);
			}
			int sindex = pid_slaveid[index];
			id[sindex].push_back(index);
		}//for
		cout<<"hello2"<<endl;
		//send random seed
		seedsize = new int[numprocs];
		randseed = new int*[numprocs];
		for(int i=1;i<numprocs;i++)
		{
			randseed[i] = new int[id[i].size()];
			seedsize[i] = id[i].size();
			for(int j=0;j<id[i].size();j++)
			{
				randseed[i][j] = id[i][j];
			}
		}
		cout<<"insert"<<endl;
		for(int i=1;i<numprocs;i++)
		{
			MPI_Send(&seedsize[i],1,MPI_INT,i,100,MPI_COMM_WORLD);
		}
		for(int i=1;i<numprocs;i++)
		{
			MPI_Send(randseed[i],seedsize[i],MPI_INT,i,99,MPI_COMM_WORLD);
		}
		cout<<"lol"<<endl;
		cout<<"hello3"<<endl;
	}
	else
	{
		MPI_Status status_a,status_b;
		int rec_size = 0;
		MPI_Recv(&rec_size,1,MPI_INT,0,100,MPI_COMM_WORLD,&status_a);
		int*slavebuf = new int[rec_size];// receive rand seed
		MPI_Recv(slavebuf,rec_size,MPI_INT,0,99,MPI_COMM_WORLD,&status_b);
		for(int i=0;i<rec_size;i++)
		{
			myarea->first_infected.insert(slavebuf[i]);
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);

	// N TIMES TRANSMIT SAME SEED
	for(int i=0;i<n;i++)
	{
		Initial_subarea();  //每个区域初始化
		MPI_Barrier(MPI_COMM_WORLD);
		for(int t=0;t<T;t++)       // transmit time from 0 to t
		{
			if(myid!=0)
			{
				myarea->just_infected.clear();
				myarea->getPrevention.clear();
	            set_iter iter = myarea->infected_now.begin();
	            set_iter iter1 = myarea->infected_now.end();
	            vector<pair<int,double> >id_age;
	            while(iter!=iter1)
	            {
	                Individual* inv_iter = myarea->Individuals.find(*iter)->second;
	                Transmit(myarea,inv_iter,t);//每个感染个体完成传播
					iter ++;
	            };
	            int_double_iter iter_id = myarea->just_infected.begin();
				int_double_iter iter_id1 = myarea->just_infected.end();
		        while(iter_id != iter_id1)//区分，也就是处理本区域的，然后存储其它区域的
		        {
		            map_iter iters = myarea->Individuals.find(iter_id->first);
		            if(iters!=myarea->Individuals.end()) // local 
		            {
		                if(iters->second->disease_state==0)
		                {
		                    myarea->infected_now.insert(iters->first);  //id
		                    iters->second->GetDisease(myarea,iter_id->second);//已经修改了diseasestate
		                }
		            }
		            else
		            {
		            	id_age.push_back(make_pair(iter_id->first,iter_id->second));
		            }
		            iter_id ++;
		        }//while,将其它的时间和id存起来
		        //处理本区域的免疫事件 vector getPrevention
		        vector<int>& getP = myarea->getPrevention;
		        int psize = getP.size();
		        int pindex = 0;
		        for(int j=0;j<psize;j++)
		        {
		        	int id_p = getP[j]; //id 
		        	map_iter iters = myarea->Individuals.find(id_p);
		        	if(iters!=myarea->Individuals.end())
		        	{
		        		if((!iters->second->is_prevention)&&iters->second->disease_state==0)
		        		{
		        			iters->second->is_prevention = true;
		        			iters->second->susceptible *= susceptible_prevention;
		        			iters->second->infectious *= infectious_prevention;
		        			iters->second->infectious_ratio *= ill_rate_prevention;
		        		}
		        	}
		        	else
		        	{
		        		prevention_id1[pindex++] = id_p;
		        	}
		        }
		        // we should send prevention_id1 to master
		        send_data_to_master(id_age,pindex);
		        slave_receive_from_master();
			}//if myid !=0
			else
			{
				master_receive_from_slave();
			}//myid==0
			MPI_Barrier(MPI_COMM_WORLD);
			//结果收集
			if(myid!=0)
			{
				int *sendbuf = new int[3];
				sendbuf[0] = myarea->dead_num;
				sendbuf[1] = myarea->recovery_num;
				sendbuf[2] = myarea->infected_now.size();
				MPI_Send(sendbuf,3,MPI_INT,0,92,MPI_COMM_WORLD);

				delete sendbuf;
			}
			else
			{
				int sum_recover = 0,sum_dead = 0,sum_infect = 0;
				int**recv_result = new int*[numprocs];
				for(int i=0;i<numprocs;i++)
				{
					recv_result[i] = new int[3];
				}
				MPI_Status status;
                for(int i=1;i<numprocs;i++)
                {
                    MPI_Recv(recv_result[i],3,MPI_INT,i,92,MPI_COMM_WORLD,&status);
					sum_recover += recv_result[i][1];
					sum_dead += recv_result[i][0];
					sum_infect += recv_result[i][2];
				}

				cout<<t<<" "<<sum_dead<<" "<<sum_recover<<" "<<sum_infect<<endl;
				ffout<<t<<" "<<sum_dead<<" "<<sum_recover<<" "<<sum_infect<<endl;
				for(int i=0;i<3;i++)
				{
					delete[] recv_result[i];
				}
				delete[] recv_result;
			}//if myid==0
			MPI_Barrier(MPI_COMM_WORLD);

		}//for t	
	} //for i
	if(myid==0)  //释放内存
	{
		cout<<"finish transmit"<<endl;
        ffout.close();
		for(int i=1;i<numprocs;i++)
		{
			delete[] masterbuf1[i];
			delete[] masterbuf2[i];
			delete[] masterbuf3[i];
			delete[] masterbuf4[i];
			delete[] masterbuf5[i];
		}
		delete[] masterbuf1;
		delete[] masterbuf2;
		delete[] masterbuf3;
		delete[] masterbuf4;
		delete[] masterbuf5;
	}
	else
	{
		delete[] infectious_id1;
		delete[] infectious_id2;
		delete[] infectious_time1;
		delete[] infectious_time2;
		delete[] prevention_id1;
		delete[] prevention_id2;
	}
	MPI_Finalize();
	return 0;
}//main
void Initial_subarea()
{
    if(myid!=0)
    {
        myarea->Reset();//感染者状态置0，当前感染置0，disease对象置NULL，diseasestate置0，表示易感者
        set_iter iter = myarea->first_infected.begin();
        set_iter iter1 = myarea->first_infected.end();
        while(iter!=iter1)
        {
            myarea->Individuals[*iter]->GetDisease(myarea,0);//初始时间为0,add to infected_now
            iter ++;
        }
    }
    
}
void send_data_to_master(vector<pair<int,double> >& v,int size)
{
    int slave_size = v.size();
    for(int i=0;i<slave_size;i++)
    {
        infectious_id1[i] = v[i].first;
        infectious_time1[i] = v[i].second;
    }
    MPI_Send(&slave_size,1,MPI_INT,0,98,MPI_COMM_WORLD);
    MPI_Send(infectious_id1,slave_size,MPI_INT,0,97,MPI_COMM_WORLD);
    MPI_Send(infectious_time1,slave_size,MPI_DOUBLE,0,96,MPI_COMM_WORLD);
    MPI_Send(&size,1,MPI_INT,0,50,MPI_COMM_WORLD);
    MPI_Send(prevention_id1,size,MPI_INT,0,49,MPI_COMM_WORLD);

}
void master_receive_from_slave()
{
    MPI_Status status,status1,status2;
    int*size = new int[numprocs];
    for(int j=1;j<numprocs;j++)
    {
        MPI_Recv(&size[j],1,MPI_INT,j,98,MPI_COMM_WORLD,&status);
        MPI_Recv(masterbuf1[j],size[j],MPI_INT,j,97,MPI_COMM_WORLD,&status1);
        MPI_Recv(masterbuf2[j],size[j],MPI_DOUBLE,j,96,MPI_COMM_WORLD,&status2);
    }
    int *sizes = new int[numprocs];
    memset(sizes,0,sizeof(int)*numprocs);
    for(int j=1;j<numprocs;j++)
    {
        for(int k=0;k<size[j];k++)
        {
            int sid = pid_slaveid[masterbuf1[j][k]];
            if(sid==0)
            {
                cout<<"fuck"<<endl<<flush;
            }
            masterbuf3[sid][sizes[sid]] = masterbuf1[j][k];
            masterbuf4[sid][sizes[sid]] = masterbuf2[j][k];
            sizes[sid] ++;
        }
    }//for
    MPI_Status status_a,status_b;
    for(int j=1;j<numprocs;j++)
    {
    	MPI_Recv(&size[j],1,MPI_INT,j,50,MPI_COMM_WORLD,&status_a);
    	MPI_Recv(masterbuf1[j],size[j],MPI_INT,j,49,MPI_COMM_WORLD,&status_b);
    }
    int *size2 = new int[numprocs];
    memset(size2,0,sizeof(int)*numprocs);
    for(int j=1;j<numprocs;j++)
    {
    	for(int k=0;k<size[j];k++)
    	{
    		int sid = pid_slaveid[masterbuf1[j][k]];
    		masterbuf5[sid][size2[sid]] = masterbuf1[j][k];
    		size2[sid] ++;
    	}
    }
    for(int j=1;j<numprocs;j++)
    {
        int ss = sizes[j];
        MPI_Send(&ss,1,MPI_INT,j,95,MPI_COMM_WORLD);
        MPI_Send(masterbuf3[j],ss,MPI_INT,j,94,MPI_COMM_WORLD);
        MPI_Send(masterbuf4[j],ss,MPI_DOUBLE,j,93,MPI_COMM_WORLD);
    }
    for(int j=1;j<numprocs;j++)
    {
    	int ss = size2[j];
    	MPI_Send(&ss,1,MPI_INT,j,48,MPI_COMM_WORLD);
        MPI_Send(masterbuf5[j],ss,MPI_INT,j,47,MPI_COMM_WORLD);
    }
    
    
}
void slave_receive_from_master()
{
    int size = 0;
    MPI_Status status1,status2,status3;
    MPI_Recv(&size,1,MPI_INT,0,95,MPI_COMM_WORLD,&status1);
    MPI_Recv(infectious_id2,size,MPI_INT,0,94,MPI_COMM_WORLD,&status2);
    MPI_Recv(infectious_time2,size,MPI_DOUBLE,0,93,MPI_COMM_WORLD,&status3);
    for(int i=0;i<size;i++)
    {
        map_iter iter = myarea->Individuals.find(infectious_id2[i]);
        if(iter == myarea->Individuals.end())
        {
            cout<<"bigerror"<<infectious_id2[i]<<" "<<myid<<endl;
            exit(0);
        }
        if(iter->second->disease_state==0)
        {
            myarea->infected_now.insert(infectious_id2[i]);//??同时感染用哪一个？
            iter->second->GetDisease(myarea,infectious_time2[i]);
        }
    }//for
    MPI_Status status_a,status_b;
    MPI_Recv(&size,1,MPI_INT,0,48,MPI_COMM_WORLD,&status_a);
    MPI_Recv(prevention_id2,size,MPI_INT,0,47,MPI_COMM_WORLD,&status_b);
    for(int i=0;i<size;i++)
    {
    	map_iter iter = myarea->Individuals.find(prevention_id2[i]);
    	if(iter == myarea->Individuals.end())
        {
            cout<<"bigerror"<<prevention_id2[i]<<" "<<myid<<endl;
            exit(0);
        }
       if((!iter->second->is_prevention)&&iter->second->disease_state==0)
		{
			iter->second->is_prevention = true;
			iter->second->susceptible *= susceptible_prevention;
			iter->second->infectious *= infectious_prevention;
			iter->second->infectious_ratio *= ill_rate_prevention;
		}
}
}
