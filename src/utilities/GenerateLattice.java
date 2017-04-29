package utilities;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import random_field.singleSpin;

public class GenerateLattice {
	
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
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+1]=new singleSpin(i*a+a/2, j*a, k*c+c/4, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+1);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+1]=new singleSpin(i*a+a/2, j*a, k*c+c/4, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+1);
            		}
        		
            		//~~~ initialize spin 2: ~~~
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+2] = new singleSpin(i*a+a/2, j*a+a/2, k*c+c/2, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+2);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+2] = new singleSpin(i*a+a/2, j*a+a/2, k*c+c/2, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+2);
            		}
        		
            		//~~~ initialize spin 3: ~~~
            		if (rnd_spin.nextDouble()<dilution){
            			int s=1;
            			if (rnd_spin.nextBoolean()) s=-1;
            			arr[i*Lx*Lz*4+j*Lz*4+k*4+3] = new singleSpin(i*a, j*a+a/2, k*c+(3*c)/4, s, rnd_spin.nextGaussian()*h, i*Lx*Lz*4+j*Lz*4+k*4+3);
            		}else{
	                	arr[i*Lx*Lz*4+j*Lz*4+k*4+3] = new singleSpin(i*a, j*a+a/2, k*c+(3*c)/4, 0, 0, i*Lx*Lz*4+j*Lz*4+k*4+3);
            		}
        		
            	}
            }
        }
        
        return arr;
    }

	public static void printArr(singleSpin[] arr, BufferedWriter out){
		for (int i=0; i<arr.length; i++){
			try{
				//if (arr[i].getX()==0){
				out.write(arr[i].getX()+","+arr[i].getY()+","+arr[i].getZ()+","+arr[i].getSpin()+","+arr[i].getH()+","+arr[i].getN());
				out.newLine();
				//}
			}
			catch (IOException e) {}
		}
	}
	
	public static void main(String[] args) {
		final double a=5.175, c=10.75;
		int Lx;	// lattice x-y size
		int Lz;	// lattice z size
        double dilution;	// starting dilution
        double h;
        int numOfConfigurations = 2;
        Lx=6;
        Lz=16;
        dilution=0.4;
        h=0;
        
        // get lattice parameters as command line arguments
        try {
        	Lx = Integer.parseInt(args[0]);
        	Lz = Integer.parseInt(args[1]);
            dilution = Double.parseDouble(args[2]);
            h = Double.parseDouble(args[3]);
            numOfConfigurations = Integer.parseInt(args[4]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
        
        for (int i=1;i<=numOfConfigurations;i++){
	        try {
	            BufferedWriter out = new BufferedWriter(new FileWriter("configurations" + File.separator + "config_"+dilution+"_"+Lx+"_"+i+".txt"));
	            
	            //print lattice sizes
	            out.write("Lx="+Lx);
	            out.newLine();
	            out.write("Lz="+Lz);
	            out.newLine();
	            
	            // print the lattice itself
	            printArr(generate_ising_lattice(Lx,Lz,a,c,dilution,h), out);
	            
	            out.close();
	        }
	        catch (IOException e) { System.out.println("bad file"); }
        }

	}

}
