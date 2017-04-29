package random_field;

public class PhaseTransCalc {

	public static void main(String[] args) {
		final double a=5.175e-10, c=10.75e-10;
		int Lx;	// lattice x-y size
		int Lz;	// lattice z size
        double dilution;	// starting dilution
        double h;
        double tau;	// EO parameter
        double gamma; // JEO aging prefactor
        int totalStepsEO;	// total steps for the EO algorithm
        int totalStepsBinder;	// total steps for binder ratio calculation
        
        // give starting values in case command-line input is empty
    	Lx=6;	//lattice x-y size
    	Lz=8;	// lattice z size
    	dilution = 0.3;
    	h = 0;
    	tau = 2;
    	gamma = 0;
    	totalStepsEO = 1;
    	totalStepsBinder = 1;
    	
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
            dilution = Double.parseDouble(args[2]);
            h = Double.parseDouble(args[3]);
            tau = Double.parseDouble(args[4]);
            gamma = Double.parseDouble(args[5]);
            totalStepsEO = Integer.parseInt(args[6]);
            totalStepsBinder = Integer.parseInt(args[7]);
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
        
        long startTime = System.currentTimeMillis();

        while (dilution>=0.280){
	        double m4=0, m2=0;
	        for (int i=0;i<totalStepsBinder;i++){
	        	singleSpin[] config = Main.main(Lx, Lz,  dilution, h, tau, gamma, totalStepsEO, a, c);
	        	double m = calcMagnetization(config);
	        	m4 += m*m*m*m;
	        	m2 += m*m;
	        }
	        m4 = m4/totalStepsBinder;	// average over the disorder
	        m2 = m2/totalStepsBinder;	// average over the disorder
	        double binder=0.5*(3-(m4/(m2*m2)));
	        System.out.println(dilution+","+binder);
	        dilution = dilution - 0.01;
        }
        
        System.out.println("run time (milliSec):"+(System.currentTimeMillis()-startTime));
        

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
					energy+=arr[i].calcDipolarInteraction(arr[j], actual_length, actual_height);
				}
			}
		}
		
		return energy/2;
		
	}
	
	public static double calcEnergy(singleSpin[] arr, double gamma){
		double energy=0;
		for(int i=1;i<arr.length;i++){
			energy = energy - (arr[i].getFitness() - arr[i].getK()*gamma) + arr[i].getAntifitness();
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
