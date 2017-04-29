package random_field;

public class trial {

	public static void main(String[] args) {
		final double a=5.175, c=10.75;
		int Lx;	// lattice x-y size
		int Lz;	// lattice z size
        double dilution;	// starting dilution
        double h;
        double tau;	// EO parameter
        double gamma; // JEO aging prefactor
        long totalStepsEO;	// total steps for the EO algorithm
        int totalStepsBinder;	// total steps for binder ratio calculation
        int fileNumber;	// configuration file number, received as parameter when running array jobs on cluster
        
        // give starting values in case command-line input is empty
    	Lx=6;	//lattice x-y size
    	Lz=8;	// lattice z size
    	dilution = 0.3;
    	h = 0;
    	tau = 2;
    	gamma = 0;
    	totalStepsEO = 1;
    	totalStepsBinder = 1;
    	fileNumber = 1;
    	
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
            dilution = Double.parseDouble(args[2]);
            h = Double.parseDouble(args[3]);
            tau = Double.parseDouble(args[4]);
            gamma = Double.parseDouble(args[5]);
            totalStepsEO = Long.parseLong(args[6]);
            totalStepsBinder = Integer.parseInt(args[7]);
            fileNumber = Integer.parseInt(args[8]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
        
        // print parameters
        System.out.println("Lx=" + Lx);
        System.out.println("Lz=" + Lz);
        System.out.println("x=" + dilution);
        System.out.println("h=" + h);
        System.out.println("tau=" + tau);
        System.out.println("gamma=" + gamma);
        System.out.println("a="+a);
        System.out.println("c="+c);
        System.out.println("total steps EO=" + totalStepsEO);
        System.out.println("total steps Binder=" + totalStepsBinder);
        System.out.println("configuration number:" + fileNumber);
        
        long startTime = System.currentTimeMillis();

        
        
        
    	//singleSpin[] config = Main.main(Lx, Lz,  dilution, h, tau, gamma, totalStepsEO, a, c, fileNumber);
        //double m = calcMagnetization(config);
        //printArr(config);
        
        int count=1;
        double[] energy_mag = Main.main(Lx, Lz,  dilution, h, tau, gamma, totalStepsEO, a, c, fileNumber);
        double minEnergy = energy_mag[0];
        double mag = energy_mag[1];
        
        for (int i=1;i<totalStepsBinder;i++){
        	energy_mag = Main.main(Lx, Lz,  dilution, h, tau, gamma, totalStepsEO, a, c, fileNumber);
        	
        	if (Math.abs(energy_mag[0]-minEnergy)<0.00001){
        		count++;
        	} else if (energy_mag[0]<minEnergy){
        		count=1;
        		minEnergy=energy_mag[0];
        		mag=energy_mag[1];
        	}
        }
        
        if (((double)count/totalStepsBinder)*100 >= 5){
        	System.out.println("m="+mag);
        }else{
        	System.out.println("ground state not found with confidence");
        }
        System.out.println((double)count/totalStepsBinder*100);
        
        //System.out.println("energy2="+calcEnergy2(config, a*Lx, c*Lz));
        
        
        System.out.println("run time (milliSec):"+(System.currentTimeMillis()-startTime));
        //System.out.println("m="+m);
        

	}
	
	public static void printArr(singleSpin[] arr){
		for (int i=1; i<arr.length; i++){
			System.out.print(arr[i].getX()+","+arr[i].getY()+","+arr[i].getZ()+","+arr[i].getSpin()+"\n");
		}
	}
	
	public static double calcEnergy2(singleSpin[] arr, double actual_length, double actual_height){
		double energy=0;
		for (int i=1;i<arr.length;i++){
			for (int j=1;j<arr.length;j++){
				if (i!=j && arr[i].getSpin()!=0 && arr[j].getSpin()!=0){
					energy+=arr[i].calcDipolarInteraction(arr[j], actual_length, actual_height) - arr[i].randomFieldEnergy();
				}
			}
		}
		
		return energy/2;
		
	}
	
	public static double calcEnergy(singleSpin[] arr, double gamma){
		double energy=0;
		for(int i=1;i<arr.length;i++){
			energy = energy - (arr[i].getFitness() - arr[i].getK()*gamma);
		}
		
		return energy/2;	// hamiltonian needs to be divided by 2
	}
	
	public static double calcMagnetization(singleSpin[] arr){
		double magnetization=0;
		int total=0;
		for (int i=1; i<arr.length; i++){
			magnetization += arr[i].getSpin();
			total+=Math.abs(arr[i].getSpin());
		}
		return magnetization/total;
	}

}
