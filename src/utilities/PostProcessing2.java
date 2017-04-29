package utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class PostProcessing2 {
	
	/**
	 * Extracts energy from a given file (file is assumed to be located in \results\)
	 * @return energy extracted from given file
	 */
	public static double[] getEnergyAndMagFromFile(int fileNumber, String fileName){
		double energy=0, m=0;
		try {
		    BufferedReader in = new BufferedReader(new FileReader("results" + File.separator + fileName + "." + fileNumber));
            String line = "";

            while ((line = in.readLine()) != null) {

                if (line.startsWith("Energy")) {
                    energy=Double.parseDouble(line.split("=")[1]);
                }
                if (line.startsWith("m=")) {
                    m=Double.parseDouble(line.split("=")[1]);
                    break;
                }
            }
		    
		    in.close();
		} catch (IOException e) {
			System.out.println("bad input file!");
		}
		return new double[]{energy, m};
	}
	public static void main(String[] args) {
		int N = 100;
		String fileName = null;
		double accuracy = 0;	// meaning perfect accuracy
        try {
        	N = Integer.parseInt(args[0]);
        	fileName = args[1];
        	accuracy = Double.parseDouble(args[2]);
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
		
        double[][] energy_mag = new double[N][2];
        double minE=0;
		for (int i=0;i<N;i++){	// file names start from 1
			energy_mag[i] = getEnergyAndMagFromFile(i+1, fileName);
			System.out.println(energy_mag[i][0]+","+energy_mag[i][1]);
			if (i==0)	minE=energy_mag[i][0];	// put first as minimum
			if (energy_mag[i][0]<minE){
				minE=energy_mag[i][0];
				System.out.println("New min found");
			}
		}
		
		int count=0;
		for (int i=0;i<N;i++){	// file names start from 1
			if (Math.abs(energy_mag[i][0]-minE)<=accuracy){
				count++;
			}
		}
		System.out.println("minimum energy:" + minE);
		System.out.println("percentage found:" + (100*(double)count/N) + "%");

	}

	

}
