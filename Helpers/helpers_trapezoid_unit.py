import numpy as np

def main():


	THETA = 0.1 # degrees

	theta_rad = THETA * np.pi / 180
	biso_per_cc = 80

	d_1 = 365 # [cm]
	h_1 =  45 # [cm]

	a_1 = 2 * d_1 * np.sin(theta_rad/2)
	b_1 = 2 * (d_1 + h_1) * np.sin(theta_rad/2)

	V_inboard = ((a_1+b_1)/2)**2 * h_1 
	N_inboard = V_inboard * biso_per_cc

	print(a_1, b_1)
	print(V_inboard, N_inboard) 



	d_2 = 830
	h_2 =  82 

	a_2 = 2 * d_2 * np.sin(theta_rad/2)
	b_2 = 2 * (d_2 + h_2) * np.sin(theta_rad/2)

	V_outboard = ((a_2+b_2)/2)**2 * h_2
	N_outboard = V_outboard * biso_per_cc

	print(a_2, b_2)
	print(V_outboard, N_outboard) 



if __name__ == '__main__':
	main()