# Photonic Ising Machine Simulation (C++)

A C++ simulation of the **photonic Ising machine equations** used in https://royalsocietypublishing.org/doi/pdf/10.1098/rsta.2021.0409.

---

## How to Compile and run 

Use `g++` to compile the code:

```bash
g++ photonic_ising_machine_simulation_one_run.cpp -o your_own_executable_name
```
Replace your_own_executable_name with the name you'd like for the compiled binary.
To run the executable: 
```bash
./your_own_executable_name N_spins N_iters alpha_min alpha_max dalpha beta_min beta_max dbeta Jij_filename h_filename energy_filename x_sol_filename
```
# Arguments
**N_spins**: Number of spins (qubits)

**N_iters**: Number of iterations

**alpha_min**, **alpha_max**, **dalpha**: Alpha range and step

**beta_min**, **beta_max**, **dbeta**: Beta range and step

**Jij_filename**: Path to file containing the Ising coupling matrix (J)

**h_filename**: Path to file containing the local field values (h)

**energy_filename**: Output filename for the computed energy

**x_sol_filename**: Output filename for the spin solution

## To do
- [ ] Add GPU-accelerated simulation of the same
- [ ] Add other source papers explaining the equations
- [ ] Add benchmarks obtained for MaxCut values using the dataset from [GSet (Stanford)](https://web.stanford.edu/~yyye/yyye/Gset/)
- [ ] Try Bayesian methods to explore parameters for better results

