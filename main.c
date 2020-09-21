/* COP 3502C Programming Assignment 3
This program is written by: Aleah Blain */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "leak_detector_c.h"

int myX;
int myY;

typedef struct coordinates{
	int x;
	int y;
	int distance;
}coordinates;

// Calcluates distance between two points
int distance(coordinates c, coordinates c2){
	return pow((c.x - c2.x),2) + pow((c.y - c2.y),2);
}

// Reads the coordinates from the input file.
coordinates *ReadData(FILE *fp, int size){
	coordinates *data = (coordinates*) malloc(size*sizeof(coordinates));

	coordinates myPoint = {myX, myY};

	for(int i = 0; i < size; i++){
		fscanf(fp, "%d %d", &data[i].x, &data[i].y);
		data[i].distance = distance(data[i], myPoint);
	}
	return data;
}

// Compares two coordinates via their distances
int compareTo(coordinates *ptrPt1, coordinates *ptrPt2){

	if(ptrPt1->distance != ptrPt2->distance)
		return ptrPt1->distance - ptrPt2->distance;

	if(ptrPt1->x != ptrPt2->x)
		return ptrPt1->x - ptrPt2->x;

	if(ptrPt1->y != ptrPt2->y)
		return ptrPt1->y - ptrPt2->y;

	return 0;
}

// Conducts a binary search within a given array for a specified coordinate
int binarySearch(coordinates c[], coordinates query, int len){
	int l = 0, r = len;    
	int mid; 
	coordinates myPoint = {myX, myY};
	int distance2;
	distance2 = distance(query, myPoint);
	query.distance = distance2;
	while (l <= r) {        
		mid = (l + r) / 2; 

		// Check if item is present at mid  
		if(compareTo(&c[mid], &query) == 0)	
			return mid; 

		// If item greater, ignore left half        
		if (compareTo(&c[mid], &query) < 0) 
			l = mid + 1;

		// If item is smaller, ignore right half       
	  else if (compareTo(&c[mid], &query) > 0)          
		  r = mid - 1;    
		}  

	// not present
	return -1;
}

// Sorts a list of coordinates via insertion method
void insertionSort(coordinates c[], int start, int end){

	int i, j;
	coordinates temp;
	for (i = start + 1; i <= end; i++)  {    
		temp = c[i];  
		for(j=i-1; j>= start; j--){  
			if(compareTo(&c[j], &temp) > 0)    
				c[j+1] = c[j];              
			else                
				break;          
		}          
	c[j+1] = temp;
	}
}

// Merges contents of an array after being sorted
void merge(coordinates c[], int start, int m, int end){
	int i, j, k;    
	int n1 = m - start + 1;    
	int n2 =  end - m;   

	/* create temp arrays */    
	coordinates *L = (coordinates*) malloc(n1*sizeof(coordinates));    
	coordinates *R = (coordinates*) malloc(n2*sizeof(coordinates));    
	
	/* Copy data to temp arrays L[] and R[] */    
	for (i = 0; i < n1; i++) {     
		L[i] = c[start + i]; 
	} 
	for (j = 0; j < n2; j++) {  
		R[j] = c[m + 1 + j];  
}
		
	/* Merge the temp arrays back into arr[l..r]*/    
	i = 0; // Initial index of first subarray    
	j = 0; // Initial index of second subarray    
	k = start; // Initial index of merged subarray    
	while (i < n1 && j < n2)    {        
		if (compareTo(&L[i], &R[j]) < 0) {            
			c[k] = L[i];             
			i++;        
		} else {            
			c[k] = R[j];            
			j++;        
		}        
		k++;   
	}   

	/* Copy the remaining elements of L[], if there are any */
	while (i < n1) {        
		c[k] = L[i]; 
		i++;        
		k++;    
	}  

	/* Copy the remaining elements of R[], if there are any */
	while (j < n2)  {        
		c[k] = R[j];       
		j++;        
		k++;    
	}
	free(L);
	free(R);
}

/* start is for left index and end is right index of the sub-array of arr to be sorted */
void mergeSort(coordinates c[], int start, int end, int t){
	if((end - start) <= t){
		insertionSort(c, start, end);
	} else {
		if (start < end) {  
			// get the mid point        
			int m = (end + start)/2; 
			// Sort first and second halves
			mergeSort(c, start, m, t); 
			mergeSort(c, m+1, end, t);
			merge(c, start, m, end);  
		}
	}	
}

// Wrapper funciton
void sort(coordinates c[], int size, int t){
		mergeSort(c, 0, size-1, t);
}

int main(void) {
	atexit(report_mem_leak);
	int size, queries, threshold;
	FILE *fp = fopen("in.txt", "r");
	FILE *fp2 = fopen("out.txt", "w");
	fscanf(fp, "%d %d %d %d %d", &myX, &myY, &size, &queries, &threshold);
	coordinates *data = ReadData(fp, size);   
	mergeSort(data, 0, size, threshold);

	// Prints sorted coordinates to out.txt
	for(int j = 0; j < size; j++)
			fprintf(fp2, "%d %d\n", data[j].x, data[j].y);

	// Searches for query within coordinates list
	for(int i = 0;  i < queries; i++){
		coordinates query;
		fscanf(fp, "%d %d", &query.x, &query.y);
		int index = binarySearch(data, query, size);
		if(index >= 0)
			fprintf(fp2, "%d %d is found at rank %d\n", query.x, query.y, index + 1);
		else
			fprintf(fp2, "%d %d not found\n", query.x, query.y);
	}
	free(data);
}