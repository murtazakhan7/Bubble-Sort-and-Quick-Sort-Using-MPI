// 2022402 task1
// sequential and paralel bublle sort

#incldue<stdio.h>
#include<time.h>
#incldue<mpi.h>

void bubblesort(int arr[], int n){
for(int i=0;i<n-1;i++){
  for(int j=0;j<n-1-i;j++{
    if(arr[j] > arr[j+1]){
      int temp = arr[j]
      arr[j] = arr[j+1];
      arr[j+1] = temp;
}
  }
}

}
int main(){
int arr[] ={1,2,23,334,0,21,4,3};
  int n = sizeof(arr)/sizeof(arr[0]);

bubblesort(arr[],n);
  for (int i=0;i<n;i++){
printf("%d ",arr[i]);
  }



return 0;
}

