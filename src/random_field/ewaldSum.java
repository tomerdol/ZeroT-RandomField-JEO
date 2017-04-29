package random_field;
import org.apache.commons.math3.special.Erf;


public class ewaldSum {
	
	public static double calcSum(singleSpin i, singleSpin j, double Lz, double Lx, double alpha, double beta){
		
		int realN_cutoff=2, reciprocalN_cutoff=7;
		double realSum=0, reciprocalSum=0;
		for (int n=0;n<=realN_cutoff;n++){
			if (n!=0 || i!=j){	// exclude self interaction
				double rz=(j.getZ()+n*Lz)-i.getZ();
				double r = Math.sqrt((j.getX()-i.getX())*(j.getX()-i.getX()) + (j.getY()-i.getY())*(j.getY()-i.getY()) + rz*rz);
				
				double B = (Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*Math.exp(-alpha*alpha*r*r))/(r*r*r);
				double C = (3*Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*(3+2*alpha*alpha*r*r)*Math.exp(-alpha*alpha*r*r))/(Math.pow(r, 5));
				
				realSum += B - rz*rz*C;
				
				// repeat with n=-n (only if n!=0)
				if (n!=0){
					rz=(j.getZ()-n*Lz)-i.getZ();
					r = Math.sqrt((j.getX()-i.getX())*(j.getX()-i.getX()) + (j.getY()-i.getY())*(j.getY()-i.getY()) + rz*rz);
					
					B = (Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*Math.exp(-alpha*alpha*r*r))/(r*r*r);
					C = (3*Erf.erfc(alpha*r) + (2*alpha*r/Math.sqrt(Math.PI))*(3+2*alpha*alpha*r*r)*Math.exp(-alpha*alpha*r*r))/(Math.pow(r, 5));
					
					realSum += B - rz*rz*C;
				}
			}
		}
		
		
		for (int mx=0;mx<=reciprocalN_cutoff/beta;mx++){
			for (int my=0;my<=reciprocalN_cutoff/beta;my++){
				for (int mz=1;mz<=reciprocalN_cutoff;mz++){
					double G2 = 4*Math.PI*Math.PI*(mx*mx*beta*beta/(Lx*Lx) + my*my*beta*beta/(Lx*Lx) + mz*mz/(Lz*Lz));
					// x positive, y positive
					reciprocalSum += (Math.pow(2*Math.PI*mz/Lz, 2)/G2)*Math.exp(-G2/(4*alpha*alpha))*2*Math.cos(2*Math.PI*mx*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*my*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*mz*(j.getZ()-i.getZ())/Lz);
					// x positive, y negative
					reciprocalSum += (Math.pow(2*Math.PI*mz/Lz, 2)/G2)*Math.exp(-G2/(4*alpha*alpha))*2*Math.cos(2*Math.PI*mx*beta*(j.getX()-i.getX())/Lx - 2*Math.PI*my*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*mz*(j.getZ()-i.getZ())/Lz);
					// x negative, y positive
					reciprocalSum += (Math.pow(2*Math.PI*mz/Lz, 2)/G2)*Math.exp(-G2/(4*alpha*alpha))*2*Math.cos(-2*Math.PI*mx*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*my*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*mz*(j.getZ()-i.getZ())/Lz);
					// x negative, y negative
					reciprocalSum += (Math.pow(2*Math.PI*mz/Lz, 2)/G2)*Math.exp(-G2/(4*alpha*alpha))*2*Math.cos(-2*Math.PI*mx*beta*(j.getX()-i.getX())/Lx - 2*Math.PI*my*beta*(j.getX()-i.getX())/Lx + 2*Math.PI*mz*(j.getZ()-i.getZ())/Lz);
				}
			}
			
		}
		
		
		
		return realSum + (4*Math.PI/(Lz*Lx*Lx))*reciprocalSum;
	}
	
	public static double realCalcSum(singleSpin i, singleSpin j, double Lz, double Lx){
		double sum = 0;
		int N_cutoff=1000;
		for (int n=0;n<=N_cutoff;n++){
			if (i!=j || n!=0){
				double rz=(j.getZ()+n*Lz)-i.getZ();
				double r = Math.sqrt((j.getX()-i.getX())*(j.getX()-i.getX()) + (j.getY()-i.getY())*(j.getY()-i.getY()) + rz*rz);
				
				sum += (r*r - 3*rz*rz)/Math.pow(r, 5);
				
				if (n!=0){
					rz=(j.getZ()-n*Lz)-i.getZ();
					r = Math.sqrt((j.getX()-i.getX())*(j.getX()-i.getX()) + (j.getY()-i.getY())*(j.getY()-i.getY()) + rz*rz);
					
					sum += (r*r - 3*rz*rz)/Math.pow(r, 5);
				}
			}
		}
		
		return sum;
	}
}
