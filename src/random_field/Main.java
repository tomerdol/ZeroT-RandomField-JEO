package random_field;

import java.util.Random;
import random_field.singleSpin;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Math;
import java.math.RoundingMode;
import java.text.DecimalFormat;


public class Main {
	
	static singleSpin[] generate_ising_lattice(int Lx, int Lz, double a, double c, double dilution, double h) {
        int i, j,k;
        // create the array that will hold the lattice. the array's cells correspond
        // to the unit cells of the LiHo{x}Y{x-1}F4.
        singleSpin[] arr = new singleSpin[4*Lx*Lx*Lz];
        Random rnd_spin = new Random();
        
        for (i = 0; i < Lx; i++)
        {
            for (j = 0; j < Lx; j++)
            {
            	for (k=0;k < Lz;k++)
            	{
            		// the spins in each unit cell are designated 0-3. see documentation for further info
            		//~~~ initialize spin 0: ~~~
            		
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+0]=new singleSpin(i*a, j*a, k*c, s,rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+0);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+0]=new singleSpin(i*a, j*a, k*c, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+0);
            		}
        		
            		//~~~ initialize spin 1: ~~~
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+1]=new singleSpin(i*a, j*a+a/2, k*c+c/4, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+1);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+1]=new singleSpin(i*a, j*a+a/2, k*c+c/4, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+1);
            		}
        		
            		//~~~ initialize spin 2: ~~~
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+2] = new singleSpin(i*a+a/2, j*a+a/2, k*c-c/2, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+2);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+2] = new singleSpin(i*a+a/2, j*a+a/2, k*c-c/2, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+2);
            		}
        		
            		//~~~ initialize spin 3: ~~~
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+3] = new singleSpin(i*a+a/2, j*a, k*c-c/4, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+3);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+3] = new singleSpin(i*a+a/2, j*a, k*c-c/4, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+3);
            		}
        		
            	}
            }
        }
        return arr;
    }
	
	static singleSpin[] receiveIsingLattice(int fileNumber, double dilution, int Lx){
		singleSpin[] arr = null;
		
		try {
		    BufferedReader in = new BufferedReader(new FileReader("configurations" + File.separator + "config_"+dilution+"_"+Lx+"_"+fileNumber+".txt"));
		    String str;
		    String[] params;
		    int i=0;
		    int Lz=0;
		    
		    if ((str = in.readLine()) != null)
	        	Lx=Integer.parseInt(str.split("=")[1]);
		    if ((str = in.readLine()) != null)		    
		    	Lz=Integer.parseInt(str.split("=")[1]);
		    
		    arr = new singleSpin[4*Lx*Lx*Lz];
		    
		    while ((str = in.readLine()) != null){
	        	params = str.split(",");
	        	arr[i]=new singleSpin(Double.parseDouble(params[0]),Double.parseDouble(params[1]), Double.parseDouble(params[2]),Integer.parseInt(params[3]),Double.parseDouble(params[4]),Integer.parseInt(params[5]));
	        	i++;
		    }
		    in.close();
		} catch (IOException e) {
			System.out.println("bad input file!");
		}
		
		return arr;
	}
	
	public static void updateAllFitness(Heap lattice, double actual_length, double actual_height, double gamma, double[][] intTable){
		singleSpin[] arr = lattice.heapArray();
		
		int i;
		for (i=1;i<arr.length;i++) {
	        int j;
	        if (arr[i].getSpin()!=0){
		    	for(j=i;j<arr.length;j++){
					if (arr[j].getSpin()!=0){
						
						
						// calculate interaction
						double currInteraction = intTable[arr[i].getN()][arr[j].getN()]*arr[i].getSpin()*arr[j].getSpin();	// dipolar and exchange interaction is kept in intTable and is then multiplied by both spin values
						arr[i].setFitness(arr[i].getFitness()  - currInteraction);
		    	    	// while we're at it, add the same interaction value to the fitness of the j spin to save us from calculating the same dipolar interaction twice
		    	    	arr[j].setFitness(arr[j].getFitness()  - currInteraction);
					}
		        }
		    	// add random field and aging parameter to fitness
		    	arr[i].setFitness(arr[i].getFitness() + arr[i].randomFieldEnergy() + gamma*arr[i].getK());
	        }
			
		}
	}
	
	public static double updateFitnessAfterFlip(Heap lattice, singleSpin i, double actual_length, double actual_height, double gamma, double[][] intTable){
		singleSpin[] arr = lattice.heapArray();
		
		i.setFitness(i.randomFieldEnergy() + gamma*i.getK());	// reset i's fitness
		double deltaEnergy=0;
    	for(int j=1;j<arr.length;j++){
			if (arr[j].getSpin()!=0){
				
				
				// get dipolar interaction from table
				double currInteraction = intTable[i.getN()][arr[j].getN()]*i.getSpin()*arr[j].getSpin();	// dipolar interaction is kept in intTable and is then multiplied by both spin values 
				i.setFitness(i.getFitness()  - currInteraction);
				
				if (i!=arr[j]){
					// while we're at it, add the same interaction value to the fitness of the j spin to save us from calculating the same dipolar interaction twice
					arr[j].setFitness(arr[j].getFitness()  - 2*currInteraction);
					deltaEnergy+=currInteraction;
				}
			}
	        
    	}
    	
    	return 2*(deltaEnergy);	//return the energy difference due to the spin flip
    	
	}
	
	public static double calcEnergy2(singleSpin[] arr, double[][] intTable, double actual_length, double actual_height){
		double energy=0;
		for (int i=0;i<arr.length;i++){
			for (int j=i;j<arr.length;j++){
				if (arr[i].getSpin()!=0 && arr[j].getSpin()!=0){
					//energy+=arr[i].calcDipolarInteraction(arr[j], actual_length, actual_height) - arr[i].randomFieldEnergy()*2;
					energy += arr[i].getSpin()*arr[j].getSpin()*intTable[arr[i].getN()][arr[j].getN()];
				}
			}
			energy += arr[i].randomFieldEnergy();
		}

		return energy;
		
	}
	
	
	// generate a random number according to the following discrete probability distribution: P(i) = 2^(-(t-1)*i)
	// n is the maximum i
	public static int generateDiscreteRandomNumber(double t, int n) {
		double[] probabilities = new double[n];
		double sum=0;
		Random rnd = new Random();
		for (int i=0;i<probabilities.length;i++){
			probabilities[i]=Math.pow(2,-(t-1)*(i+1));
			sum+=probabilities[i];
		}
		double index = rnd.nextDouble()*sum;
		double accSum=0;
		int k=0;
		while (accSum<=index){
			accSum+=probabilities[k];
			k++;
		}
		return k;
	}
	
	//calculates exchange interaction with nearest neighbors
	public static double[][] exchangeInt(double[][] intTable, singleSpin[] arr, int Lx, int Lz, double J_ex){
		for (int i=0;i<Lx;i++){	
			for (int j=0;j<Lx;j++){
				for (int k=0;k<Lz;k++){
					// notice we are only going through the 0th and 2nd atoms in the base since they participate in all exchange interactions
					
					// nearest neighbors to 0th base atom, including periodic boundary conditions
					int focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+0;
					int neighbor1=-1, neighbor2=-1, neighbor3=-1, neighbor4=-1;
					if (arr[focusSpin].getSpin()!=0){
						neighbor1=i*Lx*Lz*4+j*Lz*4+k*4+1;
				        if (i==0)
				        	neighbor2=(Lx-1)*Lx*Lz*4+j*Lz*4+k*4+1;
				        else
				        	neighbor2=(i-1)*Lx*Lz*4+j*Lz*4+k*4+1;
				        if (k==0)
				        	neighbor3=i*Lx*Lz*4+j*Lz*4+(Lz-1)*4+3;
				        else
				        	neighbor3=i*Lx*Lz*4+j*Lz*4+(k-1)*4+3;
				        
				        if (j==0 && k==0)
				        	neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(Lz-1)*4+3;
				        else if(j==0){
				        	neighbor4=i*Lx*Lz*4+(Lx-1)*Lz*4+(k-1)*4+3;
				        }
				        else if(k==0){
				        	neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(Lz-1)*4+3;
				        }else{
				        	neighbor4=i*Lx*Lz*4+(j-1)*Lz*4+(k-1)*4+3;
				        }
				        
				        // put interactions in intTable
				        if (arr[neighbor1].getSpin()!=0){
					        intTable[focusSpin][neighbor1]+=J_ex;
					        intTable[neighbor1][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor2].getSpin()!=0){
					        intTable[focusSpin][neighbor2]+=J_ex;
					        intTable[neighbor2][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor3].getSpin()!=0){
					        intTable[focusSpin][neighbor3]+=J_ex;
					        intTable[neighbor3][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor4].getSpin()!=0){
					        intTable[focusSpin][neighbor4]+=J_ex;
					        intTable[neighbor4][focusSpin]+=J_ex;
				        }
					}
			        // now nearest neighbors to 2nd base atom, including periodic boundary conditions
			        focusSpin = i*Lx*Lz*4+j*Lz*4+k*4+2;
			        if (arr[focusSpin].getSpin()!=0){
				        neighbor1=i*Lx*Lz*4+j*Lz*4+k*4+1;
				        neighbor2=i*Lx*Lz*4+j*Lz*4+k*4+3;
				        if (i==Lx-1)
				        	neighbor3=(0)*Lx*Lz*4+j*Lz*4+k*4+3;
				        else
				        	neighbor3=(i+1)*Lx*Lz*4+j*Lz*4+k*4+3;
				        if (j==Lx-1)
				        	neighbor4=i*Lx*Lz*4+(0)*Lz*4+k*4+1;
				        else
				        	neighbor4=i*Lx*Lz*4+(j+1)*Lz*4+k*4+1;

				        // put interactions in intTable
				        if (arr[neighbor1].getSpin()!=0){
					        intTable[focusSpin][neighbor1]+=J_ex;
					        intTable[neighbor1][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor2].getSpin()!=0){
					        intTable[focusSpin][neighbor2]+=J_ex;
					        intTable[neighbor2][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor3].getSpin()!=0){
					        intTable[focusSpin][neighbor3]+=J_ex;
					        intTable[neighbor3][focusSpin]+=J_ex;
				        }
				        if (arr[neighbor4].getSpin()!=0){
					        intTable[focusSpin][neighbor4]+=J_ex;
					        intTable[neighbor4][focusSpin]+=J_ex;
				        }
			        }
				}
			}
		}
		
		return intTable;
	}
	
	public static double[][] fillIntTable(double[][] intTable, singleSpin[] arr, double actual_length, double actual_height, double a, double alpha, double beta){
		double D=a*a*a*0.214;
		for (int i=0;i<arr.length;i++){
			for (int j=i;j<arr.length;j++){
				if (arr[i].getSpin()==0 || arr[j].getSpin()==0){
					intTable[i][j]=0;	// main diagonal, and missing spins are 0
				}else{
					
					double interaction = ewaldSum.calcSum(arr[i], arr[j], actual_height, actual_length, alpha, beta);
					//double interaction = ewaldSum.realCalcSum(arr[i], arr[j], actual_height, actual_length);
					
					if (i==j){	// subtract self interaction arising from the reciprocal lattice summation
						interaction -= 2*Math.pow(alpha, 3)/(3*Math.sqrt(Math.PI));
					}
					
					// ****** calculate interaction with 4 x-y plane image copies of the lattice ******
					// x+1
					double rz=arr[i].getZ()-arr[j].getZ();
					
					double r=Math.pow(arr[i].getX() + actual_length - arr[j].getX(),2) + Math.pow(arr[i].getY()-arr[j].getY(),2) + rz*rz;
					
					interaction += (r*r-3*rz*rz)/Math.pow(r, 5);
							
					// x-1
					r=Math.pow(arr[i].getX() - actual_length - arr[j].getX(),2) + Math.pow(arr[i].getY()-arr[j].getY(),2) + rz*rz;
					interaction += (r*r-3*rz*rz)/Math.pow(r, 5);
					
					// y+1
					r=Math.pow(arr[i].getX() - arr[j].getX(),2) + Math.pow(arr[i].getY() + actual_length - arr[j].getY(),2) + rz*rz;
					interaction += (r*r-3*rz*rz)/Math.pow(r, 5);
					
					// y-1
					r=Math.pow(arr[i].getX() - arr[j].getX(),2) + Math.pow(arr[i].getY() - actual_length - arr[j].getY(),2) + rz*rz;
					interaction += (r*r-3*rz*rz)/Math.pow(r, 5);
					
					// ********************************************************************************
					intTable[i][j] = D*interaction;
					intTable[j][i] = D*interaction;	// table is symmetric
					
					
					if (i==8 && j==80){
						System.out.println("interactions is: "+D*interaction);
					}
					
					
				}
			}
			
		}
		/*
		for (int i=0;i<intTable.length;i++){
			for (int j=0;j<intTable[i].length;j++){
				System.out.print(intTable[i][j]+" ");
			}
			System.out.println();
		}
		*/
		return intTable;
	}
	
	/**
	 * Runs full JEO algorithm to find ground state at zero temperature
	 * @param L lattice size
	 * @param dilution holmium atom dilution
	 * @param h random field strength
	 * @param tau EO parameter
	 * @param gamma JEO aging prefactor
	 * @param totalSteps total steps for the EO algorithm
	 * @param a lattice constant a
	 * @param c lattice constant c
	 * @return Returns the (fractional) magnetization of the system at the ground state
	 */
	public static double[] main(int Lx, int Lz, double dilution, double h, double tau, double gamma, long totalSteps, double a, double c, int fileNumber) {
		
		// final double a=5.175e-10, c=10.75e-10;
		
		/*
        int L;	// lattice size
        double dilution, h;
        double tau;	// EO parameter
        double gamma; // JEO aging prefactor
        int totalSteps;	// total steps for the EO algorithm
        */
		
        /*
        // give starting values in case command-line input is empty
    	L=6;	//lattice size
    	dilution = 0.3;
    	h = 0;
    	tau = 2;
    	gamma = 0;
    	totalSteps = 1;
    	
        try {
        	L = Integer.parseInt(args[0]);
            dilution = Double.parseDouble(args[1]);
            h = Double.parseDouble(args[2]);
            tau = Double.parseDouble(args[3]);
            gamma = Double.parseDouble(args[4]);
            totalSteps = Integer.parseInt(args[5]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
        */
        
		/*
        // print parameters
        System.out.println("L=" + L);
        System.out.println("x=" + dilution);
        System.out.println("h=" + h);
        System.out.println("tau=" + tau);
        System.out.println("gamma=" + gamma);
        System.out.println("a="+a);
        System.out.println("c="+c);
        System.out.println("total steps=" + totalSteps);
        
        */
        
        // main program
        
        // singleSpin[] bestConfig = generate_ising_lattice(Lx, Lz, a, c, dilution, h);	// the starting configuration is assumed to be best and will later be replaced
        singleSpin[] bestConfig = receiveIsingLattice(fileNumber, dilution, Lx);	// the starting configuration is assumed to be best and will later be replaced
        
        double[][] intTable = new double[4*Lx*Lx*Lz][4*Lx*Lx*Lz]; // create interaction table that holds all the dipolar interactions will be full even though it's symmetric
        intTable = fillIntTable(intTable, bestConfig, a*Lx, c*Lz, a, 1/(c*Lz), 0.5);
        intTable = exchangeInt(intTable, bestConfig, Lx, Lz, 0.12);	// add exchange interaction (J_ex is 0.12)
        
        Heap lattice = new Heap(bestConfig, tau);	// store the lattice in a heap data structure
        updateAllFitness(lattice, a*Lx, c*Lz, gamma, intTable);
        
        double minEnergy = calcEnergy(lattice.heapArray(), gamma);
        double energy = minEnergy;
        
        lattice.buildHeap();
        
        for (long step=0;step<totalSteps;step++){
        	
	        lattice.updateHeap();	// updates the heap structure with the new fitness values
        	
	        //double tempEnergy = calcEnergy2(lattice.heapArray(), Lx*a, Lz*c);
	        
	    	singleSpin flippedSpin = lattice.getRandomSpin(tau);
	    	flippedSpin.flipSpin();	// flip a random spin from the chosen level
	    	//System.out.println(step + "actual energy difference:" + (double)(calcEnergy2(lattice.heapArray(),a*Lx, c*Lz) - tempEnergy));
	    	double tempEnergy2 = updateFitnessAfterFlip(lattice, flippedSpin, a*Lx, c*Lz, gamma, intTable);
	    	//System.out.println(step + "supposed energy difference:" + tempEnergy2);
	    	energy += tempEnergy2;	// the energy difference due to the spin flip is returned
	    	//System.out.println(step + "supposed energy:" + energy);
	    	//System.out.println(step + "actual energy:" + calcEnergy2(lattice.heapArray(), Lx*a, Lz*c));
	    	
	    	if (energy < minEnergy){
	    		bestConfig = copyLattice(lattice.heapArray());
	    		minEnergy = energy;
	    	}
	    	
        }
        //System.out.println("var:"+calcSTDFitness(lattice.heapArray()));
        
        //***** print energy up to 7 decimal points to disregard floating point arithmetic errors ****
        
        //DecimalFormat df = new DecimalFormat("#.#######");
        //df.setRoundingMode(RoundingMode.HALF_DOWN);
        System.out.println("Energy="+minEnergy);
        
        // *******************************************************************************************
        //return bestConfig;
        
        printArr(bestConfig);
        //System.out.println("real energy="+calcEnergy2(bestConfig, intTable, a*Lx, c*Lz));
        return new double[]{minEnergy, calcMagnetization(bestConfig)};
    	
    	
	}
	
	public static double calcSTDFitness(singleSpin[] arr){
		double varF=0, meanF=0;
		int total =0;
		int averageK=0;
		for (int i=1; i<arr.length; i++){
			if (arr[i].getSpin()!=0){
				meanF += arr[i].getFitness();
				total++;
				averageK+=arr[i].getK();
			}
		}
		
		meanF=((double)meanF)/total;
		for (int i=1; i<arr.length; i++){
			if (arr[i].getSpin()!=0){
				varF += Math.pow(arr[i].getFitness()-meanF,2);
				
			}
		}
		//System.out.println("average k:"+((double)averageK/total));
		return Math.sqrt(((double)varF)/total);
	}
	
	public static double calcMagnetization(singleSpin[] arr){
		double magnetization=0;
		int total=0;
		for (int i=0; i<arr.length; i++){
			magnetization += arr[i].getSpin();
			total+=Math.abs(arr[i].getSpin());
		}
		return magnetization/total;
	}
	
	public static void printArr(singleSpin[] arr){
		for (int i=1; i<arr.length; i++){
			System.out.print(arr[i].getX()+","+arr[i].getY()+","+arr[i].getZ()+","+arr[i].getSpin()+"\n");
		}
	}
	

	
	// calculate the total energy (hamiltonian) by summing the fitness parameters of every spin, removing the aging parameter
	public static double calcEnergy(singleSpin[] arr, double gamma){
		double energy=0;
		for(int i=1;i<arr.length;i++){
			
			energy = energy - arr[i].getFitness() + arr[i].getK()*gamma;
			
		}
		
		return (energy/2);	// hamiltonian needs to be divided by 2
	}
	
	// copies the contents of the given iterator to an array of size N
	// used for saving a configuration
	public static singleSpin[] copyLattice(singleSpin[] arr){
		singleSpin[] newLatticeArr = new singleSpin[arr.length-1];
		for (int i=1;i<arr.length;i++){
			newLatticeArr[i-1] = new singleSpin(arr[i]);
		}
		return newLatticeArr;
	}
	
	

}
