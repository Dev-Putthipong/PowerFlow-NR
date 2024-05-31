# MATLAB Code for Power flow analysis by using Newton-Raphson method  
Created by Putthipong Niyomkitjakankul  | May 2024  

***Disclaimer : This MATLAB Code was created out of the personal interest and is intended for educational purposes only. Please note that it may contain error.***

*filename : [PowerFlowAnalysis_NR.m](PowerFlowAnalysis_NR.m)*

## Program Input

Import some test case data to MATLAB code, it includes

**BusData** : This table contains information about each bus in the power system.
*	Bus : This column identifies the number assigned to each bus. (Required)
*	Type : This column specifies the type of bus. It can be: (Slack bus, PV bus and Load/PQ bus). (Required)
*	V : This column represents the bus voltage magnitude (in per unit). Initially unknown, it can be assigned NaN or an estimated starting guess.
*	Phase : This column represents the voltage phase angle (in radians). Initially unknown, it can be assigned NaN or an estimated starting guess.
*	Pload : This column specifies the real power consumed as a load at the bus (in MW). (Required)
*	Qload : This column specifies the reactive power consumed as a load at the bus (in Mvar). (Required)
*	Pgen : This column represents the real power generation at the bus (in MW). 
*	Qgen : This column represents the reactive power generation at the bus (in Mvar). 
*	Qinjected : This column can be ignored (approximate value Q injected by shunt susceptance B) 
*	G : This column specifies the shunt conductance at the bus (Required).
*	B : This column specifies the shunt susceptance at the bus (Required).
*	Qmax : This column specifies the maximum reactive power generation limit for the generator (in Mvar) (Required for PV buses).
*	Qmin : This column specifies the minimum reactive power generation limit for the generator (in Mvar) (Required for PV buses).

**BranchData** : This table contains information about the branches (lines and transformers) connecting the buses in the power system.
*	bus_i – bus_j : These 2 columns identify the two buses (by their numbers) that a branch connects. (Required)
*	R : This column specifies the resistance of the branch (in per unit). (Required)
*	X : This column specifies the reactance of the branch (in per unit). (Required)
*	B : This column specifies line capacitive susceptance of the line (in per unit). For transformers, set this value to 0. (Required)
*	half_B : This column is B divided by 2 (Required)
*	Tap : This column specifies the tap setting of the transformer at the bus_i. For lines, set this value to 1. (Required)
*	Shift : This column specifies the phase shift angle (in degrees) for a phase. If no phase shifter is present, set this value to 0. (Required)

## Program Output
The program provides results in two ways: Workspace and Command Prompt  

### Workspace 

**Bus Information**: This table summarizes the state of each bus in the power system after convergence. It includes:
*	No_Bus : This column shows bus number
*	V : This column shows bus voltage magnitude (per unit)
*	Phase : This column shows voltage phase angle (radians)
*	P_gen : This column shows real power generation (MW)
*	Q_gen : This column shows reactive power generation (Mvar)
*	P_load : This column shows real power load (MW)
*	Q_load : This column shows reactive power load (Mvar)
*	Q_injected : This column shows reactive power injected to bus (Mvar) (due to shunt susceptance B at bus)
*	P_bus_loss : real power bus losses (MW) (due to shunt conductance G at bus)

**I_line_flow** : This table shows the current flowing through each branch (using equivalent pi model) in the network. It includes:
*	i, j : These 2 columns show index of bus number
* Ii0, Ij0 : These 2 columns show current flow out from bus i to ground and bus j to ground, respectively.
*	Iij, Iji : These 2 columns show current flow out from bus i to bus j and bus j to bus i, respectively.
*	Iij_line = Iij – Ii0 and Iji_line = Iji – Ij0
*	for phase shifter, it will show only Iij and Iji

**S_line_flow** : This table shows the complex power flowing through each branch 
*	i, j : These 2 columns show index of bus number
*	Si0, Sj0 : These 2 columns show complex power flow out from bus i to ground and bus j to ground, respectively.
*	Sij, Sji : These 2 columns show complex power flow out from bus i to bus j and bus j to bus i, respectively.
*	Sij_line = Sij – Si0 and Sji_line = Sji – Sj0
*	Sij_line+Sji_line 
*	Si0+Sj0 : total complex power flow to ground (For line, It is total line charging)
*	Total_loss : total complex power loss in the branch
*	for phase shifter, it will show only Sij, Sji and Total_loss

**Summary_line_flow** : This table provides summary of power flow of each branch in network (reduced form of S_line_flow)
*	i, j : These 2 columns show index of bus number
*	Pij, Pji : These 2 columns show real power flow out bus i to bus j and bus j to bus i, respectively.
*	Qij, Qji : These 2 columns show reactive power flow out bus i to bus j and bus j to bus i, respectively.
*	P_loss : These 2 columns show total real power loss in the branch
*	Q_loss : These 2 columns show total reactive power loss in the branch

All tables above, last row of table is total of each column.

**Y** : bus attmittance matrix   
**Psch, Qsch** : scheduled real and reactive power for calculation by Newton’s Raphson (Don’t serious to these elements, this is for who want to see procedure in detail)  
**Newton** : structure array that contains several variables of each Newton-Raphson iteration. It includes Jacobian matrix, Pcal, Qcal, Residual, Mismatch, V and Phase. The row of structure array represents round iteration. (Don’t serious to these elements, this is for who want to see procedure in detail)

### Command Prompt  
By default, the command prompt displays: Summar_line_flow and Bus Information

## Check generator limit  
If the program detects a generator exceeding its reactive power limits, it will display an alert in the command prompt. You will then be prompted to choose:
1) Fix : The program will attempt to adjust generator reactive power outputs to stay within limits and continue the analysis.
2) Don't Fix : The program will stop the analysis and display the current results.
After adjusting and repeat procedure again, it will check bus voltage magnitude if it is suitable.  If it is not suitable, it will display an alert in the command prompt
and prompt you for a similar fix/don't fix decision. Once the program successfully completes the analysis without any violations, it will display "The power flow analysis is completed"
and you can see the results.

## Noted
* ***I downloaded power system test case data from the University of Washington Power Systems Test Case Archive ([UW Power Systems Test Case Archive](https://labs.ece.uw.edu/pstca/)). I have formatted the data for use in my code verification process. This is for educational purposes only. Please note that the formatted data which I provided may contain errors during the formatting process.***
*	You can adjust S_base, tolerance and initial guess for NaN in code editor
*	To import your data, please keep the original names of your data including both table and column names.
*	For easier organization, I recommend using Excel to prepare your data before importing.
