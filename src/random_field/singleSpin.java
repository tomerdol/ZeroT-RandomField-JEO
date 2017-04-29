package random_field;


public class singleSpin implements Comparable<singleSpin>{
	private double x,y,z;	// coordinates
	private int s;			// spin
	private double h;		// local random field
	private double fitness;	// EO fitness
	private int k;			// aging parameter
	private int n;			// ordinal number
	
	// 1st constructor, which sets the position to the default (0,0,0)
	public singleSpin(){
		this(0, 0, 0, 1, 0, 0);
	}
	
	//2nd constructor, which sets the position to a given (x0,y0,z0)
	public singleSpin(double x0, double y0, double z0, int s0, double h0, int n0){
		x=x0;
		y=y0;
		z=z0;
		h=h0;
		n=n0;
		if (s0==1 || s0==-1 || s0==0){
			s=s0;
		} else{
			throw new IllegalArgumentException("Invalid spin! spin should be 1 or -1 or 0");
		}
		
		fitness = 0;	// fitness will be determined only after the lattice is complete
		k = 0;			// the aging parameter records how many times a spin has been selected by the EO algorithm. obviously starts at 0
		
	}
	
	
	public singleSpin(singleSpin other){
		this(other.x, other.y, other.z, other.s, other.h, other.n);
	}
	
	// getters to return the required values
	public double getX(){
		return this.x;
	}
	public double getY(){
		return this.y;
	}
	public double getZ(){
		return this.z;
	}
	public int getSpin(){
		return this.s;
	}
	public double getH(){
		return this.h;
	}
	public double getFitness(){
		return this.fitness;
	}
	public int getK(){
		return this.k;
	}
	public int getN(){
		return this.n;
	}
	
	//setters
	public void setX(double x){
		this.x=x;
	}
	public void setY(double y){
		this.y=y;
	}
	public void setZ(double z){
		this.z=z;
	}
	public void setSpin(int s0){
		if (s0==1 || s0==-1 || s0==0){
			this.s=s0;
		} else{
			throw new IllegalArgumentException("Invalid spin! spin should be 1 or -1");
		}
	}
	public void setH(double h){
		this.h=h;
	}
	public void setFitness(double fitness){
		this.fitness=fitness;
	}
	public void setK(int k){
		this.k=k;
	}
	public void setN(int n){
		this.n=n;
	}
	//toString prints the spin & position of the spin in the following form: 'Spin:-1 Position:(x,y,z)'
	public String toString(){
		return "Spin:"+s+" Position:("+this.x+","+this.y+","+this.z+")" + " Local Random Field:"+h+" Fitness:"+fitness;
	}
	
	//returns the distance between this point and another
	public double distance(singleSpin another){
		return Math.sqrt(Math.pow(this.x-another.x, 2)+Math.pow(this.y-another.y, 2)+Math.pow(this.z-another.z, 2));
	}
	
	//distance function for periodic boundary conditions
	public double distance(singleSpin another, double xL, double zL){
		double dx = this.x - another.x;
		double dy = this.y - another.y;
		double dz = this.z - another.z;
		// implementation of the periodic boundary conditions for each dimension separately
		if (dx>xL*0.5)
			dx=dx-xL;
		if (dx<=-xL*0.5)
			dx=dx+xL;
		if (dy>xL*0.5)
			dy=dy-xL;
		if (dy<=-xL*0.5)
			dy=dy+xL;
		if (dz>zL*0.5)
			dz=dz-zL;
		if (dz<=-zL*0.5)
			dz=dz+zL;

		return Math.sqrt(dx*dx+dy*dy+dz*dz);
	}
	
	//flips the spin
	public void flipSpin(){
		s=s*(-1);	// flip the spin
		k=k+1; // add 1 to the aging parameter
	}
	
	// returns the energy from the spin's interaction with the local random field
	public double randomFieldEnergy(){
		return this.s*this.h;
	}
	
	//compares 2 spins according to their fitness
	public int compareTo(singleSpin other){
		if (this.fitness - other.fitness < 0)
			return -1;
		else if (this.fitness - other.fitness > 0)
			return 1;
		else
			return 0;
		
	}
	
	public double calcDipolarInteraction(singleSpin o, double actual_length, double actual_height){
		if (this==o)	return 0;
		double r=this.distance(o, actual_length, actual_height);
		double rz=this.z-o.z;
		//double rx=i.getX()-j.getX();
		
		// implementation of periodic boundary conditions
		if (rz>actual_height*0.5)
			rz=rz-actual_height;
		if (rz<=-actual_height*0.5)
			rz=rz+actual_height;
		
		return (this.s*o.s)*(r*r-3*rz*rz)/Math.pow(r, 5);
		
	}


	
	
	
	
}
