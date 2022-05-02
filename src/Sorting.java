import java.util.Random;
import Plotter.Plotter;


public class Sorting{

	final static int SELECT_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int COUNTING_VS_QUICK_LENGTH = 16;
	final static int BUBBLE_VS_MERGE_LENGTH = 12;
	final static int MERGE_VS_QUICK_SORTED_LENGTH =11;
	final static double T = 600.0;

	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 *
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 *
	 * @param arr - the array to be sorted
	 */
	public static void quickSort(double[] arr, int firstIndex, int lastIndex){
		if (firstIndex < lastIndex){
			int pivotIndex = partition(arr, firstIndex, lastIndex);
			quickSort(arr, firstIndex, pivotIndex-1);
			quickSort(arr, pivotIndex+1, lastIndex);
		}
	}

	private static int partition(double[] arr, int firstIndex, int lastIndex) {
		double pivot = arr[lastIndex];
		int i = (firstIndex-1);
		for (int j = firstIndex; j < lastIndex; j++) {
			if (arr[j] <= pivot) {
				i++;
				swap(arr, i, j);
			}
		}

		swap(arr, i+1, lastIndex);
		return i+1;
	}

	private static void swap(double[] arr, int firstIndex, int lastIndex){
		double temporaryValue = arr[firstIndex];
		arr[firstIndex] = arr[lastIndex];
		arr[lastIndex] = temporaryValue;
	}

	/**
	 * Given an array arr and an index i returns the the i'th order statstics in arr.
	 * In other words, it returns the element with rank i in the array arr.
	 *
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 *
	 * Should run in average complexity of O(n), and worst case complexity of O(n^2)
	 *
	 **/
	public static double QuickSelect(double[] arr, int firstIndex, int lastIndex, int rank) {
		if (firstIndex == lastIndex) return arr[firstIndex];
		int pivotIndex = partition(arr, firstIndex, lastIndex);
		int pivotRank = pivotIndex - firstIndex + 1;
		if (rank == pivotRank) return arr[pivotIndex];
		else if (rank < pivotRank){
			return QuickSelect(arr, firstIndex, pivotIndex - 1, rank);
		}
		else {
			return QuickSelect(arr, pivotIndex + 1, lastIndex, rank - pivotRank);
		}
	}

	/**
	 * Sorts a given array using the merge sort algorithm.
	 * 
	 * Should run in complexity O(nlog(n)) in the worst case.
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr, int firstIndex, int lastIndex){
		if (firstIndex < lastIndex) {
			int medianIndex = firstIndex + (lastIndex - firstIndex) / 2;
			mergeSort(arr, firstIndex, medianIndex);
			mergeSort(arr, medianIndex+1, lastIndex);
			merge(arr, firstIndex, medianIndex, lastIndex);
		}
	}

	private static void merge(double[] arr, int firstIndex, int medianIndex, int lastIndex) {
		int leftSubarraySize = medianIndex - firstIndex + 1;
		int rightSubarraySize = lastIndex - medianIndex;
		double[] leftSubarray = new double[leftSubarraySize];
		double[] rightSubarray = new double[rightSubarraySize];

		System.arraycopy(arr, firstIndex, leftSubarray, 0, leftSubarraySize);
		for (int j = 0; j < rightSubarraySize; j++) {
			if (medianIndex + 1 + j >= arr.length)  return;
			else {
				rightSubarray[j] = arr[medianIndex + 1 + j];
			}
		}

		addElementsBySize(arr, leftSubarray, rightSubarray, rightSubarraySize, leftSubarraySize, firstIndex);
	}

	private static void addElementsBySize(double[] arr, double[] leftSubarray, double[] rightSubarray,int rightSubarraySize, int leftSubarraySize, int firstIndex) {
		int indexOfLeftSubarray = 0;
		int indexOfRightSubarray = 0;
		int indexOfMergedArray = firstIndex;
		while ((indexOfLeftSubarray < leftSubarraySize) && (indexOfRightSubarray < rightSubarraySize)) {
			if (leftSubarray[indexOfLeftSubarray] <= rightSubarray[indexOfRightSubarray]) {
				arr[indexOfMergedArray] = leftSubarray[indexOfLeftSubarray];
				indexOfLeftSubarray++;
			} else {
				arr[indexOfMergedArray] = rightSubarray[indexOfRightSubarray];
				indexOfRightSubarray++;
			}
			indexOfMergedArray++;
		}
		addRestSubarrayElements(arr, leftSubarray, leftSubarraySize, indexOfMergedArray, indexOfLeftSubarray);
		addRestSubarrayElements(arr, rightSubarray, rightSubarraySize, indexOfMergedArray, indexOfRightSubarray);
	}

	private static void addRestSubarrayElements(double[] arr, double[] subarray, int subarraySize, int indexOfArray, int indexOfSubarray){
		while (indexOfSubarray < subarraySize) {
			arr[indexOfArray] = subarray[indexOfSubarray];
			indexOfSubarray++;
			indexOfArray++;
		}
	}

	/**
	 * Sorts a given array using bubble sort.
	 * 
	 * The algorithm should run in complexity O(n^2).
	 * 
	 * @param arr - the array to be sorted
	 */
	public static void bubbleSort(double[] arr){
		boolean stop = false;
		for (int i = 0; i < arr.length - 2; i++){
			if (stop) return;
			boolean swap = false;
			for (int j = 1; j < arr.length - i; j++){
				if (arr[j-1] > arr[j]) {
					swap = true;
					double temp = arr[j-1];
					arr[j-1] = arr[j];
					arr[j] = temp;
				}
			}

			if (!swap) stop = true;
		}
	}
	/**
	 * Sorts a given array, using the counting sort algorithm.
	 * You may assume that all elements in the array are between 0 and k (not including k).
	 * 
	 * Should run in complexity O(n + k) in the worst case.
	 * 
	 * @param arr - an array with possitive integers
	 * @param k - an upper bound for the values of all elements in the array.
	 */
	public static void countingSort(int[] arr,int[] sortedArr, int k){
		int[] arrRanks = new int[k];
		for(int i = 0; i < arr.length; i++){
			arrRanks[arr[i]] += 1;
		}
		for(int j = 1; j < arrRanks.length; j++){
			arrRanks[j] = arrRanks[j] + arrRanks[j-1];
		}
		for(int r = arr.length - 1; r >= 0; r--){
			sortedArr[arrRanks[arr[r]]-1] = arr[r];
			arrRanks[arr[r]] -= 1;
		}
	}

	public static void main(String[] args) {
		countingVsQuick();
		mergeVsQuick();
		mergeVsQuickOnSortedArray();
		mergeVsBubble();
		QuickSelectVsQuickSort();
	}

	private static void countingVsQuick() {
		double[] quickTimes = new double[COUNTING_VS_QUICK_LENGTH];
		double[] countingTimes = new double[COUNTING_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < COUNTING_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumCounting = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				int[] b = new int[size];
				for (int j = 0; j < a.length; j++) {
					b[j] = r.nextInt(size);
					a[j] = b[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a, 0, a.length-1);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				int[] sortedB = new int[b.length];
				countingSort(b, sortedB, size);
				endTime = System.currentTimeMillis();
				sumCounting += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			countingTimes[i] = sumCounting/T;
		}
		Plotter.plot("Counting sort on arrays with elements < n", countingTimes, "Quick sort on arrays with elements < n", quickTimes);
		
	}
	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a,0, a.length-1);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b, 0, b.length-1);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
	}
	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a,0, a.length-1);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b, 0, b.length - 1);
				endTime = System.currentTimeMillis();
				sumMerge  += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}
	/**
	 * Compares merge sort and bubble sort on random arrays
	 */
	public static void mergeVsBubble(){
		double[] mergeTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		double[] bubbleTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < BUBBLE_VS_MERGE_LENGTH; i++) {
			long sumMerge = 0;
			long sumBubble = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				mergeSort(a, 0, a.length);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				bubbleSort(b);
				endTime = System.currentTimeMillis();
				sumBubble += endTime - startTime;
			}
			mergeTimes[i] = sumMerge/T;
			bubbleTimes[i] = sumBubble/T;
		}
		Plotter.plot("merge sort on random array", mergeTimes, "bubble sort on random array", bubbleTimes);
	}
	/**
	 * Compares the quick select algorithm with a random rank, and the quick sort algorithm.
	 */
	public static void QuickSelectVsQuickSort(){
		double[] QsortTimes = new double[SELECT_VS_QUICK_LENGTH];
		double[] QselectTimes = new double[SELECT_VS_QUICK_LENGTH];
		Random r = new Random();
		long startTime, endTime;
		for (int i = 0; i < SELECT_VS_QUICK_LENGTH; i++) {
			long sumQsort = 0;
			long sumQselect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a,0, a.length-1);
				endTime = System.currentTimeMillis();
				sumQsort += endTime - startTime;
				startTime = System.currentTimeMillis();
				QuickSelect(b,0, b.length - 1, r.nextInt(size)+1);
				endTime = System.currentTimeMillis();
				sumQselect += endTime - startTime;
			}
			QsortTimes[i] = sumQsort/T;
			QselectTimes[i] = sumQselect/T;
		}
		Plotter.plot("quick sort with an arbitrary pivot", QsortTimes, "quick select with an arbitrary pivot, and a random rank", QselectTimes);
	}
	

}
