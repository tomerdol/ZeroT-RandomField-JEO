package utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class PostProcessing {
	
	
	/**
	 * Extracts magnetization from a given file (file is assumed to be located in \results\)
	 * @return magnetization extracted from given file
	 */
	public static double getMagnetizationFromFile(int fileNumber, String fileName){
		double m=0;
		try {
		    BufferedReader in = new BufferedReader(new FileReader("results" + File.separator + fileName + "." + fileNumber));
		    String str, lastLine=null;
		    
		    while ((str=in.readLine()) != null){
                if (str.startsWith("m=")) {
                    m=Double.parseDouble(str.split("=")[1]);
                    break;
                }
		    }
		    
		    in.close();
		} catch (IOException e) {
			System.out.println("bad input file!");
		}
		return m;
	}
	
	public static void main(String[] args) {
		int Nsa = 100;
		String fileName = null;
        try {
        	Nsa = Integer.parseInt(args[0]);
        	fileName = args[1];
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("ArrayIndexOutOfBoundsException caught");
        }
        catch (NumberFormatException e){}
		
        double m2=0, m4=0;
        
		for (int i=1;i<=Nsa;i++){	// file names start from 1
			double m = getMagnetizationFromFile(i, fileName);
			m2 += m*m;
			m4 += m*m*m*m;
		}
		
		m2=m2/Nsa;
		m4=m4/Nsa;
		double binder=0.5*(3-(m4/(m2*m2)));
		System.out.println(binder);

	}

}
