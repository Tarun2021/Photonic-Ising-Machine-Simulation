#include <bits/stdc++.h>
using namespace std;
void load_Jij(const string& filename,std::vector<std::vector<double>>& Jij, int& m,int& n) //to load the coupling matrix
{
    ifstream file(filename.c_str());
    if(!file)
    {
        cerr<<"Error: could not open file!"<<endl;
        exit(1);
    }
    file>>m>>n;    
	Jij.resize(m, std::vector<double>(n));	
    for(int i=0;i<m;i++)
    {
		for(int j=0;j<n;j++)
		{
        file>>Jij[i][j];		 
		}
    }    
    file.close();
}
void load_h(const string& filename,std::vector<double>& h, int& m) // to load the bias values (strength of external magnaetic field) for each spin
{
	ifstream file(filename.c_str());
	if(!file)
	{
		cerr<<"Error: could not open file!"<<endl;
		exit(1);
	}
	file>>m;
	h.resize(m); //needed to resize the vector to hold m elements
	for(int i=0;i<m;i++)
	{
		file>>h[i];
	}		
	file.close();
}
void initialize_array(double* arr, int N_spins) {
    srand(time(0));  // Seed for random number generation
    for (int j = 0; j < N_spins; j++) {
        arr[j] = -0.1 + static_cast<double>(rand()) / RAND_MAX * 0.2; 
    }
}
int main(int argc, char** argv) {
	int N_spins=atoi(argv[1]);
	int Niters=atoi(argv[2]);
	double alpha_min=atof(argv[3]);
	double alpha_max=atof(argv[4]);
	double daplha=atof(argv[5]);
	double beta_min=atof(argv[6]);
	double beta_max=atof(argv[7]);
	double dbeta=atof(argv[8]);	
	string Jij_filename=argv[9];
	string h_filename=argv[10];
	string energy_filename=argv[11];
	string x_sol_filename=argv[12];
	srand(time(0));	
	vector<vector<double>> Jij;
	vector<double> h;
	load_Jij(Jij_filename, Jij, N_spins, N_spins); // Load Jij matrix from file
    load_h(h_filename, h, N_spins); // Load h matrix from file
	int nzj=0;
	int nzh=0;
	double max_Jij=0.0;
	double max_h=0.0;
	for(int i=0;i<N_spins;i++)
	{
		for(int j=0;j<N_spins;j++)
		{
			if(Jij[i][j]!=0.0)
			{
				nzj++;
				if(abs(Jij[i][j])>max_Jij)
				{
					max_Jij=abs(Jij[i][j]);
				}
			}
		}
	}
	nzj/=2; //since Jij is symmetric, we count only half of the matrix
	for(int i=0;i<N_spins;i++)
	{
		if(h[i]!=0.0)
		{
			nzh++;
		}
		if(abs(h[i])>max_h)
		{
			max_h=abs(h[i]);
		}
	}
	int steps_alpha=round((alpha_max-alpha_min)/daplha);	
	int steps_beta=round((beta_max-beta_min)/dbeta);	
	double x_sol_alpha[steps_alpha][N_spins];
	for(int i=0;i<steps_alpha;i++)
	{
		initialize_array(x_sol_alpha[i], N_spins);
	}
	double alpha_t_j[steps_alpha];
	double beta_t_j[steps_beta];	
	for (int i=0;i<steps_alpha;i++)
	{
	    alpha_t_j[i]=alpha_min+i*daplha;
	}
	for (int i=0;i<steps_beta;i++)
	{
	    beta_t_j[i]=beta_min+i*dbeta;
	}
	double energy[steps_alpha][Niters];
	for(int i=0;i<steps_alpha;i++)
	{
		for(int j=0;j<Niters;j++)
		{
			energy[i][j]=0.0;
		}
	}
	
	for(int step=0;step<steps_alpha;step++) // for various values of alpha, keeping beta as 0 
	{
        double alpha=alpha_t_j[step];
		//find time taken to execute
		auto start = chrono::high_resolution_clock::now();		 
		double beta=N_spins*(alpha-1)/(max_Jij*nzj+max_h*nzh); // beta is 0 only for verification of the bifurcation; else beta is obtained from the user or from alpha
		cout<<"Beta value "<<beta<<endl;
		double feedback_signal[N_spins];
		for (int j = 0; j < N_spins; j++) {
			feedback_signal[j] = 0.0;  // Initialize to zero to avoid garbage values
		}
		double fpga_noise[N_spins];
		for (int j = 0; j < N_spins; j++) {
			fpga_noise[j] = 0.0;  // Initialize to zero to avoid garbage values
		}
		int fpga_noise_size=sizeof(fpga_noise)/sizeof(fpga_noise[0]);
		double* new_fpga_noise=gaussian_noise_added(fpga_noise,fpga_noise_size);		
		for (int i=0;i<Niters;i++)
		{
			for (int j=0;j<N_spins;j++)
		    {			
		        double dot_product=0.0;
		        for(int k=0;k<N_spins;k++)
		        {
					dot_product+=-(Jij[k][j]*x_sol_alpha[step][k]+h[k]);
		        }				
		        feedback_signal[j]=alpha*x_sol_alpha[step][j]+beta*dot_product; //feedback							
                x_sol_alpha[step][j]=-0.5+pow(cos(feedback_signal[j]-(M_PI/4)),2);//system dynamics	
						
		    }
			for(int ii=0;ii<N_spins;ii++)
			{
				for(int ij=0;ij<N_spins;ij++)
				{
					energy[step][i]+=(Jij[ii][ij]*(x_sol_alpha[step][ii]/abs(x_sol_alpha[step][ii]))*(x_sol_alpha[step][ij]/abs(x_sol_alpha[step][ij]))); //energy calculation
				}
			}			
			energy[step][i]=-(energy[step][i]);					    
		}			
		delete new_fpga_noise;
	}
	ofstream file;
	file.open(energy_filename.c_str());
	for(int i=0;i<steps_alpha;i++)
	{
		for(int j=0;j<Niters;j++)
		{
			file<<energy[i][j]<<" ";
		}
		file<<endl;
	}
	file.close();	
	file.open(x_sol_filename.c_str());
	for(int i=0;i<steps_alpha;i++)
	{
		for(int j=0;j<N_spins;j++)
		{
			file<<x_sol_alpha[i][j]<<" ";
			cout<<x_sol_alpha[i][j]/abs(x_sol_alpha[i][j])<<" ";
		}
		file<<endl;
		cout<<endl;
	}
	file.close();
	return 0;
 
}