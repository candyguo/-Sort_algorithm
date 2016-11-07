#include <iostream>
using namespace std;

int max(int *arr,int len)
{
	int max=0;
	for(int i=0;i<len;i++)
	{
		if(arr[i]>max)
			max=arr[i];
	}
	return max;
}

void print_arr(int *arr,int len)
{
	for(int i=0;i<len;i++)
		cout<<arr[i]<<"  "; 
}
//冒泡排序，从后往前，每次把最小的放在最前面
void Bubble_sort(int *arr,int len)
{
	for(int i=1;i<len;i++)
	{
		for(int j=len-1;j>0;j--)
		{
			if(arr[j]<arr[j-1])
			{
				int temp=arr[j];
				arr[j]=arr[j-1];
				arr[j-1]=temp;
			}
		}
	}
}
//选择排序，每次排序选择最小的放入有序区
void Select_sort(int *arr,int len)
{
	int minIndex=0;
	for(int i=1;i<len;i++)
	{
		minIndex=i-1;//每次排序的最小索引
		for(int j=i;j<len;j++)
		{
			if(arr[j]<arr[minIndex])
			{
				minIndex=j;
			}
		}
		if(minIndex!=i-1)
		{
			int temp=arr[i-1];
			arr[i-1]=arr[minIndex];
			arr[minIndex]=temp;
		}
	}
}

//快速排序是冒泡的改进，前后指针移动，基准比较
void Quick_sort(int *arr,int left,int right)
{
	if(left>right)
	 	return;

	// int left=first;
    int last=right;//保留初始端值迭代
	int std=arr[left];
	while(left<right)
	{
		while(left<right&&arr[right]>=std)
		{
			--right;
		}
		arr[left]=arr[right];

		while(left<right&&arr[left]<=std)
		{
			++left;
		}
		arr[right]=arr[left];
	}
	arr[right]=std;//左右指针相遇，把基准值移动到相遇的地方
	Quick_sort(arr,0,left-1);
	Quick_sort(arr,left+1,last);
}

//计数排序 时间复杂度o(n),但占用了空间
void Count_sort(int *arr,int len)
{
	int maxnum=max(arr,len);
	int *result=new int[len];
	int *Count = new int[maxnum + 1];
	for(int i=0;i<len;i++)
	{
		Count[arr[i]]++;//数字出现次数统计
	}
	//次数累加
	for(int j=1;j<=maxnum;j++)
	{
		Count[j]+=Count[j-1];
	}
	for(int k=len-1;k>=0;--k)
	{
		result[Count[arr[k]]-1] = arr[k];
		Count[arr[k]]--;
	}
	 for(int i=0;i<len;i++)
	 {
	 	arr[i]=result[i];
	 }
}

//堆排序，调整，建立，排序
void HeapAdjust(int *a,int i,int size)  //调整堆 
{
    int lchild=2*i;       //i的左孩子节点序号 
    int rchild=2*i+1;     //i的右孩子节点序号 
    int max=i;            //临时变量 
    if(i<=size/2)          //如果i不是叶节点就不用进行调整 
    {
        if(lchild<size&&a[lchild]>a[max])
        {
            max=lchild;
        }    
        if(rchild<size&&a[rchild]>a[max])
        {
            max=rchild;
        }
        if(max!=i)
        {
            swap(a[i],a[max]);
            HeapAdjust(a,max,size);    //避免调整之后以max为父节点的子树不是堆 
        }
    }        
}

void BuildHeap(int *a,int size)    //建立堆 
{
    int i;
    for(i=size-1;i>=0;i--)    //非叶节点最大序号值为size/2 
    {
        HeapAdjust(a,i,size);    
    }    
} 

void HeapSort(int *a,int size)    //堆排序 
{
    int i;
    BuildHeap(a,size);
    for(i=size-1;i>=0;i--)
    {
        swap(a[0],a[i]);        //交换堆顶和最后一个元素，即每次将剩余元素中的最大者放到最后面 
        HeapAdjust(a,0,i);      //重新调整堆顶节点成为大顶堆，<size,交换的最后一个数不参与比较
    }
} 

int main()
{
	int arr[]={5,3,4,8,6,9,7,10};
	int  *result=new int[10];
	print_arr(arr,sizeof(arr)/sizeof(arr[0]));
	HeapSort(arr,sizeof(arr)/sizeof(arr[0])); 
	cout<<endl;
	print_arr(arr,sizeof(arr)/sizeof(arr[0]));
}
