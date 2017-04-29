package random_field;

import java.util.Iterator;;

public class ArrayIterator<AnyType> implements Iterator<AnyType>{
	private int current;
	private AnyType[] arr;
	
	public ArrayIterator(AnyType[] a){
		arr=a;
		current=1;
	}
	
	public boolean hasNext(){
		return current < arr.length;
	}
	
	public AnyType next(){
		return arr[current++];
	}
	
	
	
	
}
