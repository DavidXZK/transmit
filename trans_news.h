#ifndef  TRANS_NEW_H
#define  TRANS_NEW_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

using namespace std;
using namespace std::tr1;
const int MAX_SLAVE = 500000;
const int SEED_NUM = 2000;
typedef unordered_map<int,int>::iterator id_age_iter;
class Individual;
int n,T,cutage;
double LAMDA,infectious_rate,death_rate;
double mean_latency,mean_infectious_pre,mean_infectious_sym,mean_infectious_asym;
double Age_infectious_rate,Age_susceptible_rate;
double Infectious_rate_pre,Infectious_rate_sym,Infectious_rate_asym;
double school_intensity,home_intensity,work_intensity;
double community_intensity,friends_intensity,commute_intensity;
double treat_infectious,susceptible_prevention,infectious_prevention;
double ill_rate_prevention,ill_rate_treat,infectious_treat;
double get_treat_rate,get_immune_rate;
double school_decrease,work_decrease;
void Parameter(ifstream& fin)
{
    fin>>n>>T>>cutage;
    fin>>LAMDA>>infectious_rate>>death_rate;
    fin>>mean_latency>>mean_infectious_pre>>mean_infectious_sym>>mean_infectious_asym;
    fin>>Age_infectious_rate>>Age_susceptible_rate;
    fin>>Infectious_rate_pre>>Infectious_rate_sym>>Infectious_rate_asym;
    fin>>school_intensity>>home_intensity>>work_intensity>>community_intensity;
    fin>>friends_intensity>>commute_intensity;
    fin>>treat_infectious>>susceptible_prevention>>infectious_prevention;
    fin>>ill_rate_prevention>>ill_rate_treat>>infectious_treat;
    fin>>get_treat_rate>>get_immune_rate;
	fin>>school_decrease>>work_decrease;
}

class Disease
{
	public:
		double time_infected;
		double latency;     //潜伏期
		double infect_pre_time;  //前驱期时间
		double infect_sym_time;  //症状期时间
		double infect_asym_time; //非症状期时间
		bool sym_or_not,is_pre_cal,is_sym_cal,is_asym_cal;//是否计算过了
		Disease(double time,double latency,double infect_pre_time,double infect_sym_time,double infect_asym_time,bool sym_or_not)
		:time_infected(time),latency(latency),infect_pre_time(infect_pre_time),infect_sym_time(infect_sym_time),infect_asym_time(infect_asym_time),sym_or_not(sym_or_not),is_pre_cal(false)
		{
			is_pre_cal = false;
			is_sym_cal = false;
			is_asym_cal = false;
		}
};

class SubArea
{
	public:
		int area_id,treat_drug_num,prevent_drug_num,dead_num,recovery_num;
		unordered_map<int,Individual*>Individuals;
		unordered_map<int,double>just_infected;
		unordered_set<int> infected_now,first_infected; //即时感染者，初始感染者
		vector<int> getPrevention;//获得免疫接种的人
		SubArea(int id):area_id(id),treat_drug_num(0),prevent_drug_num(0),dead_num(0),recovery_num(0){}
		void Reset();
};

class Individual
{
	public:
		int pid,age,slaveid,disease_state;
		double susceptible,infectious,infectious_ratio;
		bool is_treat,is_prevention;//治疗以及预防
		Disease*disease;
		unordered_map<int,int> home,school,work,friends,community,commute;  //id-age
		unordered_map<int,double> temp_infectious;  //临时感染?
		Individual(int pid,int age,int slaveid):pid(pid),age(age),slaveid(slaveid),infectious_ratio(infectious_rate)
		{
			is_treat = false;
			is_prevention = false;
			susceptible = 1.0;
			infectious =1.0;
			disease_state = 0;
		}
		void GetDisease(SubArea* apt,double time);
		void FreeDisease(SubArea* apt);
		void Dead(SubArea* apt);
		void Reset()
		{

			disease_state = 0;
			susceptible = 1.0;
			infectious = 1.0;
			is_treat = false;
			is_prevention = false;
			delete disease;
			disease = NULL;
		}
};
typedef unordered_map<int,double>::iterator int_double_iter;
typedef unordered_map<int,Individual*>::iterator map_iter;
typedef unordered_set<int>::iterator set_iter;
void readNet(string filenames,int flag,SubArea*apt);
void CalInfect(SubArea* apt,Individual* inv,double re_infectious_rate,double current_time,double period);
void Transmit(SubArea* apt,Individual* inv,double current_time);
int getImmune(SubArea* apt,Individual* inv);
void readNet(string filenames,int flag,SubArea*apt)
{
    string filename = filenames;
	ifstream fin(filename.c_str());
	int pid,degree,destid,destage;
	string line;
	while(fin>>pid>>degree)
	{
		//fin>>pid>>degree;
        unordered_map<int,Individual*>::iterator it = apt->Individuals.find(pid);
		if(it!=apt->Individuals.end())
		{
			for(int i =0;i<degree;i++)
			{
				fin>>destid>>destage;
                switch(flag)
                {
                    case 1:
                        it->second->home[destid] = destage;
                        break;
                    case 2:
                        it->second->school[destid] = destage;
                        break;
                    case 3:
                        it->second->work[destid] = destage;
                        break;
                    case 4:
                        it->second->friends[destid] = destage;
                        break;
                    case 5:
                        it->second->community[destid] = destage;
                        break;
                    case 6:
                        it->second->commute[destid] = destage;
                        break;
                    default:
                        cout<<"readNet default fault"<<endl;
                        break;
                }//switch	
			}
		}
		else
		{
			getline(fin,line);
		}
	}//while
	fin.close();
}
void SubArea::Reset()
{
	treat_drug_num = 0;
	prevent_drug_num = 0;
	dead_num = 0;
	recovery_num = 0;
	just_infected.clear();
	infected_now.clear();
	getPrevention.clear();
	map_iter iter1 = Individuals.begin();
	map_iter iter2 = Individuals.end();
	while(iter1 !=iter2)
	{
		if(iter1->second->disease_state != 0)
		{
			iter1->second->Reset();
		}
		iter1 ++;
	}
}
void Individual::GetDisease(SubArea* apt,double time)
{
	double random1 = (double)rand()/RAND_MAX;
	double random2 = (double)rand()/RAND_MAX;
	double latency = (-mean_latency)*log(1-random2);
	if(random1<infectious_ratio)// get disease
	{
		random1 = (double)rand()/RAND_MAX;
		double infect_pre = (-mean_infectious_pre)*log(1-random1);
		random2 = (double)rand()/RAND_MAX;
		double infect_sym = (-mean_infectious_sym)*log(1-random2);
		disease = new Disease(time,latency,infect_pre,infect_sym,0,true);
		disease_state = 1;
		apt->infected_now.insert(pid);
	}
	else
	{
		random1 = (double)rand()/RAND_MAX;//无症状感染期
		double infect_asy = (-mean_infectious_asym)*log(1-random1);
		disease = new Disease(time,latency,0,0,infect_asy,false);
		disease_state = 1;
		apt->infected_now.insert(pid);
	}

}
void Individual::FreeDisease(SubArea* apt)
{
	disease_state = 2;
	apt->recovery_num ++;
	apt->infected_now.erase(pid);
	delete disease;
	disease = NULL;
}
void Individual::Dead(SubArea* apt)
{
	disease_state = 3;
	apt->dead_num ++;
	apt->infected_now.erase(pid);
	delete disease;
	disease = NULL;
}
void getTreat(Individual*inv)
{
	double rate = (double)rand()/(double)RAND_MAX;
	if(rate <= get_treat_rate)
	{
		inv->is_treat = true;
		inv->infectious *= infectious_treat;
	}
}
int getImmune(SubArea* apt,Individual*inv)
{
	double rate = (double)rand()/(double)RAND_MAX;
	if(rate<=get_immune_rate)   //社会网络中的易感者接受免疫
	{
		unordered_map<int,int>::iterator iter = inv->home.begin();
		while(iter !=inv->home.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		iter = inv->school.begin();
		while(iter !=inv->school.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		iter = inv->friends.begin();
		while(iter !=inv->friends.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		iter = inv->work.begin();
		while(iter !=inv->work.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		iter = inv->community.begin();
		while(iter !=inv->community.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		iter = inv->commute.begin();
		while(iter !=inv->commute.end())
		{
			apt->getPrevention.push_back(iter->first);
			iter ++;
		}
		return 1;
	}//if
	else
	{
		return 0;
	}
}
void CalInfect(SubArea* apt,Individual* inv,double re_infectious_rate,double current_time,double period)
{
	double ageir = ((inv->age <= cutage)?Age_infectious_rate:1.0)*inv->infectious;
	id_age_iter iter1 = inv->home.begin();
	id_age_iter iter2 = inv->home.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*home_intensity*ageir;
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
	iter1 = inv->school.begin();
	iter2 = inv->school.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*school_intensity*ageir*(1-school_decrease);
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
	iter1 = inv->work.begin();
	iter2 = inv->work.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*work_intensity*ageir*(1-work_decrease);
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
	iter1 = inv->friends.begin();
	iter2 = inv->friends.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*friends_intensity*ageir;
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
	iter1 = inv->community.begin();
	iter2 = inv->community.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*community_intensity*ageir;
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
	iter1 = inv->commute.begin();
	iter2 = inv->commute.end();
	while(iter1 !=iter2)
	{
		int age = iter1->second;
		double lamda = ((age>=cutage)?1.0:Age_susceptible_rate)*
		LAMDA*re_infectious_rate*commute_intensity*ageir;
		double random = (double)rand()/RAND_MAX;
		double infect_time = (-(double)1/lamda)*log(1-random);
		if(infect_time <= period)
		{
			apt->just_infected[iter1->first] = (double)current_time + infect_time;
		}
		iter1 ++;
	}
}
void Transmit(SubArea* apt,Individual* inv,double current_time)
{
	double latency = inv->disease->latency;
	double time_infected = inv->disease->time_infected;
	double time_pre = inv->disease->infect_pre_time;
	double time_sym = inv->disease->infect_sym_time;
	double time_asym = inv->disease->infect_asym_time;
	if(inv->disease->sym_or_not)//显式感染
	{
		double sym_time = current_time - time_infected;
		if((!inv->disease->is_pre_cal)&&(latency<=sym_time))// pre
		{
			double period = inv->disease->infect_pre_time;
			double re_infectious_rate = Infectious_rate_pre;
			CalInfect(apt,inv,re_infectious_rate,current_time,period);
			inv->disease->is_pre_cal = true;
		}
		if((!inv->disease->is_sym_cal)&&(sym_time>latency+time_pre))//sym
		{
		//	if(sym_time>=latency+time_pre+time_sym)
		//	{
		//		cout<<"alert!"<<endl;
		//		return;
		//	}
			double period = inv->disease->infect_sym_time;    //症状感染期长度
			double re_infectious_rate = Infectious_rate_sym;
			if(getImmune(apt,inv)==1){
				re_infectious_rate *= susceptible_prevention; 
			}
			getTreat(inv);
			CalInfect(apt,inv,re_infectious_rate,current_time,period);
			inv->disease->is_sym_cal = true;
		//	return;
		}
		if((sym_time>=latency+time_pre+time_sym)&&inv->disease->is_pre_cal&&inv->disease->is_sym_cal)//康复或者死亡
		{
			double random = (double)rand()/RAND_MAX;
			double dr = death_rate;
			if(random<=dr)
			{
				inv->Dead(apt);
			}
			else
			{
				inv->FreeDisease(apt);
			}
		}
	}
	else  //无症状感染
	{
		double asym_time = current_time - time_infected;
		if(!inv->disease->is_asym_cal&&(asym_time>=latency))
		{
			 double period = inv->disease->infect_asym_time;
			 double re_infectious_rate = Infectious_rate_asym;
			 CalInfect(apt,inv,re_infectious_rate,current_time,period);
			 inv->disease->is_asym_cal = true;
		}
		if(inv->disease->is_asym_cal&&(asym_time>=(latency+time_asym)))//结束，康复
        {
            inv->FreeDisease(apt);
        }
	}
}

#endif
