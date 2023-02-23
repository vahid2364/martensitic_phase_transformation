# martensitic_phase_transformation
Phase field modeling of the martensitic phase transformation

<img src="icon.png" align="right" />

# Cahn-Hilliard-Equation-Solver-FDM-2D # [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)]

## filename = CH.f90

	Author : Vahid Attari
	Created: 30 Feb. 2016
	Modified: ....
	Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
	Acknowledgements:  Based on Cahn-Hilliard 1965 paper
	
## Purpose:

	- Phase Field Modeling with dynamic coupling to thermodyanmic and kinetic databases to self consistantly model the Spinodal Composition Phenomenon
	- Dynamic coupling is not provided ...
   
## General Algorithm function:

	1. Retrieve parameter data from file "parameters.dat"
	2. Assess thermodynamics of the associated system ()
	3. Reads initial phase distribution from "phase.dat" file (Not provided)
	4. Calculate Phase Evolution with time integration
	5. Nucleate Phases (Not provided)

	-  Resolve boundary conditions 
	 	-- Periodic boundaries in all directions

	-  Solve differential equations via 9-stencil finite difference
	-  Update phase information and concentration data


!! Compilation instructions: 
	>> make file is not provided
	>> - Manual: >>  ifort -o a.out CH.f90

!! Execution: >> ./a.out 
                                
	!!------------------------------------------------------------------------------------
	!!------------------------------------------------------------------------------------
	!!------------------------------------------------------------------------------------
	!!====================================================================================

	!   This code simulates the early stage of Spinodal Decomposition...

	!   Cahn- Hilliard solver 	
	!   with periodic boundary conditions.
	!   Nonlinear term: f(u) = u - u**3

<table>
  <tr>
    <td> 
<img src="https://user-images.githubusercontent.com/11892854/118386432-78937500-b5e5-11eb-9c48-dc04c4be50b4.jpeg" alt="microstructures_000001" width="250" height="250">
<img src="https://user-images.githubusercontent.com/11892854/118386435-792c0b80-b5e5-11eb-84fd-5f993fc2c2c2.jpeg" alt="microstructures_000002" width="250" height="250">
<img src="https://user-images.githubusercontent.com/11892854/118386436-7a5d3880-b5e5-11eb-915a-dd687dc01aaf.jpeg" alt="microstructures_000004" width="250" height="250">
<img src="https://user-images.githubusercontent.com/11892854/118386439-7cbf9280-b5e5-11eb-9cee-3cc17a7f0ab2.jpeg" alt="microstructures_000007" width="250" height="250">	    
<img src="https://user-images.githubusercontent.com/11892854/118386440-7d582900-b5e5-11eb-89a2-bc0252ba8135.jpeg" alt="microstructures_000009" width="250" height="250">	    
<img src="https://user-images.githubusercontent.com/11892854/118386441-7d582900-b5e5-11eb-98bb-9520f8dfe864.jpeg" alt="microstructures_000011" width="250" height="250">	    
<img src="https://user-images.githubusercontent.com/11892854/118386442-7df0bf80-b5e5-11eb-8f20-d5262b2c1b3d.jpeg" alt="microstructures_000014" width="250" height="250">	    
<img src="https://user-images.githubusercontent.com/11892854/118386443-7df0bf80-b5e5-11eb-82a5-eba0a00faeae.jpeg" alt="microstructures_000016" width="250" height="250">	    
<img src="https://user-images.githubusercontent.com/11892854/118386444-7e895600-b5e5-11eb-9bc1-b40dfb0002cd.jpeg" alt="microstructures_000021" width="250" height="250">	    
	  </td>
   <tr>
</table>
