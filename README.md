<img src="icon.png" align="right" />

# Martensitic_phase_transformation # [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)]
Phase field modeling of the martensitic phase transformation

## filename = Allen_cahn_2D.f90

	Author : Vahid Attari
	Created: 30 Feb. 2022
	Modified: ....
	Arroyave Research Group, Department of Materials Science & Engineering, Texas A&M University
	Acknowledgements:  Based on 
	
## Purpose:

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


## Compilation instructions: 
	>> see make file
	>> - Manual: >>  ifort -qmkl $MKLROOT/mkl_dfti.f90 -o pfm.out file.f90

!! Execution: >> ./a.out 
                                
	!!------------------------------------------------------------------------------------
	!!------------------------------------------------------------------------------------
	!!------------------------------------------------------------------------------------
	!!====================================================================================

	!   This code simulates the early stage of transition ...

	!   Two Allen- Hilliard solver 	
	!   with periodic boundary conditions.
	!   Nonlinear term: f(u) = 

## Watch the video:

<table>
  <tr>
    <td> 
 \url(https://github.com/vahid2364/martensitic_phase_transformation/blob/main/animation/animation.mp4)
	  </td>
   <tr>
</table>
