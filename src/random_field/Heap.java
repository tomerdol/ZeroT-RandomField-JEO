package random_field;

import java.util.Random;
import java.lang.Math;

/**
 * @author Tomer
 *
 */
public class Heap {
	private static final int CAPACITY = 2;

	private int size; // Number of elements in heap
	private singleSpin[] heap; // The heap array
	private int height;	// height of the heap
	private double[] probabilities;	// probability distribution for choosing random spins
	
	public Heap() {
		size = 0;
		height = 0;
		heap = new singleSpin[CAPACITY];
	}

	/**
	 * Construct the binary heap given an array of items.
	 */
	public Heap(singleSpin[] array, double tau) {
		
		size = 0;
		for (int i=0;i<array.length;i++){
			if (array[i].getSpin()!=0)	size++;
		}
		
		heap = new singleSpin[size + 1];
		
		int heapIndex=1;	// we do not use the index 0
		for (int i=0;i<array.length;i++){
			if (array[i].getSpin()!=0){
				heap[heapIndex]=array[i];
				heapIndex++;
			}
		}
		
		// calculate heap height
		int i = 1;
		while (Math.pow(2, i) <= size) {
			i++;
		}
		height=i;
		
		// populate probability array
		// the array has a cell for each level of the heap
		// each cell contains a number proportional to the probability of choosing that level
		probabilities = new double[height];
		double sum=0;
		for (i=0;i<probabilities.length;i++){
			probabilities[i]=Math.pow(2,-(tau-1)*(i+1));
			sum+=probabilities[i];
		}
		// normalize probability:
		for (i=0;i<probabilities.length;i++){
			probabilities[i]=probabilities[i]/sum;
		}
		
		// to finish, order the heap for the first time:
		buildHeap();
	}
	
	/**
	 * returns the size of the heap
	 */
	public int getSize(){
		return size;
	}
	
	/**
	 * runs at O(size)
	 */
	public void buildHeap() {
		for (int k = size / 2; k > 0; k--) {
			percolatingDown(k);
		}
	}
	
	/**
	 * performs a partial update similar to the one described in S. Boettcher, Extremal optimization for Sherrington-Kirkpatrick spin glasses, Eur. Phys. J. B 46, 501–505 (2005)
	 */
	public void updateHeap() {
		for (int k=1; k <= size/2; k++){
			singleSpin tmp = heap[k];
			if (k*2 <= size){
				int child=k*2;
				if (child != size && heap[child].compareTo(heap[child + 1]) > 0)
					child++;
				if (tmp.compareTo(heap[child]) > 0){
					heap[k] = heap[child];
					heap[child]=tmp;
				}
			}
		}
	}

	private void percolatingDown(int k) {
		singleSpin tmp = heap[k];
		int child;

		for (; 2 * k <= size; k = child) {
			child = 2 * k;

			if (child != size && heap[child].compareTo(heap[child + 1]) > 0)
				child++;

			if (tmp.compareTo(heap[child]) > 0)
				heap[k] = heap[child];
			else
				break;
		}
		heap[k] = tmp;
	}

	/**
	 * Sorts a given array of items.
	 */
	public void heapSort(singleSpin[] array) {
		size = array.length;
		heap = (singleSpin[]) new singleSpin[size + 1];
		System.arraycopy(array, 0, heap, 1, size);
		buildHeap();

		for (int i = size; i > 0; i--) {
			singleSpin tmp = (singleSpin) heap[i]; // move top item to the end of the heap
										// array
			heap[i] = heap[1];
			heap[1] = tmp;
			size--;
			percolatingDown(1);
		}
		for (int k = 0; k < heap.length - 1; k++)
			array[k] = (singleSpin) heap[heap.length - 1 - k];
	}

	/**
	 * Deletes the top item
	 */
	public singleSpin deleteMin() throws RuntimeException {
		if (size == 0)
			throw new RuntimeException();
		singleSpin min = (singleSpin) heap[1];
		heap[1] = heap[size--];
		percolatingDown(1);
		
		// update heap height
		int i = 1;
		while (Math.pow(2, i) <= size) {
			i++;
		}
		height=i;
		
		return min;
	}

	/**
	 * Inserts a new item
	 */
	public void insert(singleSpin x) {
		if (size == heap.length - 1)
			doubleSize();

		// Insert a new item to the end of the array
		int pos = ++size;

		// Percolate up
		for (; pos > 1 && x.compareTo((singleSpin) heap[pos / 2]) < 0; pos = pos / 2)
			heap[pos] = heap[pos / 2];

		heap[pos] = x;
		
		// update heap height
		int i = 1;
		while (Math.pow(2, i) <= size) {
			i++;
		}
		height=i;
	}

	private void doubleSize() {
		singleSpin[] old = (singleSpin[]) heap;
		heap = (singleSpin[]) new singleSpin[heap.length * 2];
		System.arraycopy(old, 1, heap, 1, size);
	}

	public String toString() {
		String out = "";
		for (int k = 1; k <= size; k++)
			out += heap[k] + " ";
		return out;
	}
	
	/**
	 * Returns a random spin from the given level
	 * @param level - the level from which to return a random spin
	 */
	public singleSpin getRandomSpin(double tau) {
		// generate a random number according to the probability in probabilities[] array
		Random rnd = new Random();
		
		double index = rnd.nextDouble();
		double accSum=0;
		int level=0;	// will be used in a sec as the level from which to choose a spin
		while (accSum<=index){
			accSum+=probabilities[level];
			level++;
		}
		
		// now get a random spin from the determined level
		int upperBound = (int) Math.pow(2, level); // upper bound, exclusive
		int lowerBound = (int) Math.pow(2, level - 1); // lower bound, inclusive
		if (upperBound - 1 > size)
			upperBound = size + 1;
		
		return heap[lowerBound + rnd.nextInt(upperBound - lowerBound)];
	}

	public int getHeight() {
		return height;
	}
	
	/*
	public Iterator iterator() {
		return new ArrayIterator((singleSpin[]) heap);
	}
	*/
	public singleSpin[] heapArray(){
		return heap;
	}

}
