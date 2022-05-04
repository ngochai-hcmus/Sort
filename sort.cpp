#include<iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <stdlib.h>
#include <string>

using namespace std;

int a[1000000], temp[1000000];
int Count[500001];

//SELECTION SORT
//https://nguyenvanhieu.vn/thuat-toan-sap-xep-selection-sort/

//RUNTIME
double SelectionSort(int arr[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    int i, j, min_idx;
    for (i = 0; i < n - 1; i++)
    {
        min_idx = i;
        for (j = i + 1; j < n; j++)
        {
            if (arr[j] < arr[min_idx])
            {
                min_idx = j;
            }
        }
        swap(arr[min_idx], arr[i]);
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void SelectionSort_CP(int arr[], int n, unsigned long long int& count_compar)
{
    count_compar = 0;
    int i = 0;
    while (++count_compar && i < n - 1)
    {
        int j = i + 1;
        while (++count_compar && j < n)
        {
            int  min_idx = i;
            if (++count_compar && arr[j] < arr[min_idx])
            {
                min_idx = j;
            }
            swap(arr[min_idx], arr[i]);
            j++;
        }
        i++;
    }
}


//INSERTION SORT
//https://www.geeksforgeeks.org/insertion-sort/

//RUNTIME
double InsertionSort(int arr[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    for (int i = 1; i < n; i++)
    {
        int value = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > value)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = value;
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void InsertionSort_CP(int arr[], int n, unsigned long long int& count_compar)
{
    count_compar = 0;
    int i = 1;
    while (++count_compar && i < n)
    {
        int value = arr[i];
        int j = i - 1;
        while (++count_compar && j >= 0 && arr[j] > value)
        {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = value;
        i++;
    }
}


//BUBBLE SORT
//https://www.geeksforgeeks.org/bubble-sort/

//RUNTIME
double BubbleSort(int arr[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    for (int i = 0; i < n - 1; ++i)
    {
        bool swapped = false;
        for (int j = 0; j < n - i - 1; ++j)
        {
            if (arr[j] > arr[j + 1])
            {
                swap(arr[j], arr[j + 1]);
                swapped = true;
            }
        }

        if (swapped == false)
            break;
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void BubbleSort_CP(int arr[], int n, unsigned long long int& count_compar)
{
    count_compar = 0;
    int i = 0;
    while (++count_compar && i < n - 1)
    {
        bool swapped = false;
        int j = 0;
        while (++count_compar && j < n - i)
        {
            if (++count_compar && arr[j] > arr[j + 1])
            {
                swap(arr[j], arr[j + 1]);
                swapped = true;
            }
            j++;
        }
        if (swapped == false)
            break;
        i++;
    }
}


//HEAP SORT

//RUNTIME
void Heapify(int arr[], int n, int i)
{
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest])
        largest = left;

    if (right < n && arr[right] > arr[largest])
        largest = right;

    if (largest != i) {
        swap(arr[i], arr[largest]);
        Heapify(arr, n, largest);
    }
}

double HeapSort(int arr[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    for (int i = n / 2 - 1; i >= 0; i--)
        Heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--)
    {
        swap(arr[0], arr[i]);
        Heapify(arr, i, 0);
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void Heapify_CP(int arr[], int n, int i, unsigned long long int& count_compar)
{
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left > n) count_compar++;
    else {
        if (left < n && arr[left] > arr[largest])
        {
            largest = left;
            count_compar += 2;
        }
    }
    if (right > n) count_compar++;
    else {
        if (right < n && arr[right] > arr[largest])
        {
            count_compar += 2;
            largest = right;
        }
    }
    if (largest != i) {
        count_compar++;
        swap(arr[i], arr[largest]);
        Heapify_CP(arr, n, largest, count_compar);
    }
}

void HeapSort_CP(int arr[], int n, unsigned long long int& count_compar)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        count_compar++;
        Heapify_CP(arr, n, i, count_compar);
    }

    for (int i = n - 1; i > 0; i--)
    {
        count_compar++;
        swap(arr[0], arr[i]);
        Heapify_CP(arr, i, 0, count_compar);
    }
}

//MERGE SORT

//RUNTIME
void Merge(int arr[], int  l, int  m, int  r)
{
    int i, j, k;
    int   n1 = m - l + 1;
    int n2 = r - m;
    int* L = new int[n1];
    int* R = new int[n2];

    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i] <= R[j])
        {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}

double MergeSort(int arr[], int l, int r)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    if (l < r)
    {
        int m = l + (r - l) / 2;
        MergeSort(arr, l, m);
        MergeSort(arr, m + 1, r);
        Merge(arr, l, m, r);
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void Merge_CP(int arr[], int l, int m, int r, unsigned long long int& count_compar)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    int* L = new int[n1];
    int* R = new int[n2];

    for (i = 0; i < n1; i++)
    {
        L[i] = arr[l + i];
        count_compar++;
    }
    for (j = 0; j < n2; j++)
    {
        R[j] = arr[m + 1 + j];
        count_compar++;
    }
    i = 0;
    j = 0;
    k = l;
    if (i > n1) count_compar++;
    else {
        while (i < n1 && j < n2)
        {
            count_compar += 2;
            if (L[i] <= R[j])
            {
                count_compar++;
                arr[k] = L[i];
                i++;
            }
            else {
                arr[k] = R[j];
                j++;
            }
            k++;
        }
        while (i < n1)
        {
            count_compar++;
            arr[k] = L[i];
            i++;
            k++;
        }
        while (j < n2)
        {
            count_compar++;
            arr[k] = R[j];
            j++;
            k++;
        }
    }
}

void MergeSort_CP(int arr[], int l, int r, unsigned long long int& count_compar)
{
    if (l < r)
    {
        count_compar++;
        int m = l + (r - l) / 2;
        MergeSort_CP(arr, l, m, count_compar);
        MergeSort_CP(arr, m + 1, r, count_compar);
        Merge_CP(arr, l, m, r, count_compar);
    }
}


//QUICK SORT

//RUNTIME
double QuickSort(int arr[], int low, int high) {
    clock_t timeStart, timeEnd;
    timeStart = clock();
    int i, j, x;
    x = arr[(low + high) / 2];
    i = low, j = high;
    do {
        while (arr[i] < x) i++;
        while (arr[j] > x) j--;
        if (i <= j) {
            swap(arr[i], arr[j]);
            i++;
            j--;
        }
    } while (i <= j);
    if (low < j) QuickSort(arr, low, j);
    if (i < high) QuickSort(arr, i, high);
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void QuickSort_CP(int arr[], int low, int high, unsigned long long int& count_compar) {
    int i, j, x;
    x = arr[(low + high) / 2];
    i = low, j = high;
    do {
        count_compar++;
        while (arr[i] < x) {
            i++;
            count_compar++;
        }
        while (arr[j] > x) {
            j--;
            count_compar++;
        }
        if (i <= j) {
            count_compar++;
            swap(arr[i], arr[j]);
            i++;
            j--;
        }
    } while (i <= j);
    if (low < j) {
        count_compar++;
        QuickSort_CP(arr, low, j, count_compar);
    }
    if (i < high) {
        count_compar++;
        QuickSort_CP(arr, i, high, count_compar);
    }
}


//RADIX SORT
//https://www.geeksforgeeks.org/radix-sort/
//RUNTIME
int getBit(int a, int k)
{
    return ((a >> k) & 1);
}

double RadixSort(int a[], int x, int y, int k)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    if (k < 0 || x >= y)
    {
        timeEnd = clock();
        double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
        return timeRun;
    }
    int i, j;
    for (i = x, j = y; i < j;)
    {
        while (i < y && getBit(a[i], k) == 0) ++i;
        while (j >= x && getBit(a[j], k) == 1) --j;
        if (i < j) swap(a[i], a[j]);
    }
    RadixSort(a, x, j, k - 1);
    RadixSort(a, i, y, k - 1);
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
int getBit_CP(int a, int k, unsigned long long int& count_compar)
{
    return ((a >> k) & 1);
}

void RadixSort_CP(int a[], int x, int y, int k, unsigned long long int& count_compar)
{
    if (++count_compar && k < 0 || x >= y)
    {
        return;
    }
    int i, j;
    for (i = x, j = y; ++count_compar && i < j;)
    {
        while (++count_compar && i < y && getBit_CP(a[i], k, count_compar) == 0) ++i;
        while (++count_compar && j >= x && getBit_CP(a[j], k, count_compar) == 1) --j;
        if (++count_compar && i < j) swap(a[i], a[j]);
    }
    RadixSort_CP(a, x, j, k - 1, count_compar);
    RadixSort_CP(a, i, y, k - 1, count_compar);
}


//SHAKER SORT
//https://www.stdio.vn/giai-thuat-lap-trinh/bubble-sort-va-shaker-sort-01Si3U

//RUNTIME
double ShakerSort(int a[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    int Left = 0;
    int Right = n - 1;
    int k = 0;
    while (Left < Right)
    {
        for (int i = Left; i < Right; i++)
        {
            if (a[i] > a[i + 1])
            {
                swap(a[i], a[i + 1]);
                k = i;
            }
        }
        Right = k;
        for (int i = Right; i > Left; i--)
        {
            if (a[i] < a[i - 1])
            {
                swap(a[i], a[i - 1]);
                k = i;
            }
        }
        Left = k;
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void ShakerSort_CP(int a[], int n, unsigned long long int& count_compar)
{
    int Left = 0;
    int Right = n - 1;
    int k = 0;
    while (++count_compar&& Left < Right)
    {
        for (int i = Left; ++count_compar && i < Right; i++)
        {
            if (++count_compar && a[i] > a[i + 1])
            {
                swap(a[i], a[i + 1]);
                k = i;
            }
        }
        Right = k;
        for (int i = Right; ++count_compar && i > Left; i--)
        {
            if (++count_compar&& a[i] < a[i - 1])
            {
                swap(a[i], a[i - 1]);
                k = i;
            }
        }
        Left = k;
    }
}


//SHELL SORT
//https://www.softwaretestinghelp.com/shell-sort/

//RUNTIME
double ShellSort(int arr[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    for (int gap = n / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i += 1) {
            int temp = arr[i];
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];
            arr[j] = temp;
        }
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void ShellSort_CP(int arr[], int n, unsigned long long int& count_compar) {
    count_compar = 0;
    int gap = n / 2;
    while (++count_compar && gap > 0) {
        int i = gap;
        while (++count_compar && i < n) {
            int temp = arr[i];
            int j = i;
            while ((++count_compar && j >= gap) && (++count_compar && arr[j - gap] > temp)) {
                arr[j] = arr[j - gap];
                j -= gap;
            }
            arr[j] = temp;
            i += 1;
        }
        gap /= 2;
    }
}


//COUNTING SORT
//https://www.geeksforgeeks.org/counting-sort/
//RUNTIME
double CountingSort(int a[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();
    for (int i = 0; i <= 100000; ++i)
        Count[i] = 0;
    for (int i = 0; i < n; ++i)
        ++Count[a[i]];
    int temp = 0;
    for (int i = 0; i < n; ++i)
    {
        while (Count[temp] <= 0)
            ++temp;
        a[i] = temp;
        --Count[temp];
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void CountingSort_CP(int a[], int n, unsigned long long int& count_compar)
{
    count_compar = 0;
    //int Count[100001]; // a[i] <= 100000
    for (int i = 0; ++count_compar && i <= 100000; ++i)
        Count[i] = 0;
    for (int i = 0; ++count_compar && i < n; ++i)
        ++Count[a[i]];
    int temp = 0;
    for (int i = 0; ++count_compar && i < n; ++i)
    {
        while (++count_compar && Count[temp] <= 0)
            ++temp;
        a[i] = temp;
        --Count[temp];
    }
}


//FLASH SORT
//https://www.youtube.com/watch?v=CAaDJJUszvE

//RUNTIME
double FlashSort(int a[], int n)
{
    clock_t timeStart, timeEnd;
    timeStart = clock();

    int m = int(0.45 * n);// 0.45 là tu cong thưc cua flash sort ma ra
    for (int i = 0; i < n; ++i)
    {
        if (a[0] < a[i]) swap(a[0], a[i]); // max=a[0]
        if (a[n - 1] > a[i]) swap(a[n - 1], a[i]); //min = a[n-1]
    }

    int mn = a[n - 1];
    double tempValue = double(m - 1) / (a[0] - a[n - 1]);
    vector<int> buckets(m, 0);
    for (int i = 0, k; i < n; ++i)
    {
        k = int(tempValue * (a[i] - mn)); //Cong thưc cua flash sort: (m-1) * (a[i]-min) / (max - min)
        ++buckets[k];
    }

    for (int i = 1; i < m; ++i)
        buckets[i] += buckets[i - 1];

    for (int i = 0, j = 0, k = m - 1; i < n - 1;)
    {
        if (k < 0) break;
        while (j >= buckets[k]) //ban dau a[0] = max nen se khong gap truong hop nay
        {
            ++j;
            k = int(tempValue * (a[j] - mn));
        }

        while (j != buckets[k])
        {
            k = int(tempValue * (a[j] - mn));
            swap(a[j], a[buckets[k] - 1]);
            --buckets[k];
            ++i;
        }
    }

    for (int i = 1, j; i < n; ++i) //insertion sort
    {
        j = i;
        while (a[j] < a[j - 1])
        {
            swap(a[j], a[j - 1]);
            --j;
        }
    }
    timeEnd = clock();
    double timeRun = double(timeEnd - timeStart) / double(CLOCKS_PER_SEC);
    return timeRun;
}

//COUNT COMPARE
void FlashSort_CP(int a[], int n, unsigned long long int& count_compar)
{
    count_compar = 0;

    int m = int(0.45 * n);// 0.45 là tu cong thưc cua flash sort ma ra
    for (int i = 0; ++count_compar && i < n; ++i)
    {
        if (++count_compar && a[0] < a[i]) swap(a[0], a[i]); // max=a[0]
        if (++count_compar && a[n - 1] > a[i]) swap(a[n - 1], a[i]); //min = a[n-1]
    }

    int mn = a[n - 1];
    double tempValue = double(m - 1) / (a[0] - a[n - 1]);
    vector<int> buckets(m, 0);
    for (int i = 0, k; ++count_compar && i < n; ++i)
    {
        k = int(tempValue * (a[i] - mn)); //Cong thưc cua flash sort: (m-1) * (a[i]-min) / (max - min)
        ++buckets[k];
    }

    for (int i = 1; ++count_compar && i < m; ++i)
        buckets[i] += buckets[i - 1];

    for (int i = 0, j = 0, k = m - 1; ++count_compar && i < n - 1;)
    {
        if (++count_compar && k < 0) break;
        while (++count_compar && j >= buckets[k]) //ban dau a[0] = max nen se khong gap truong hop nay
        {
            ++j;
            k = int(tempValue * (a[j] - mn));
        }

        while (++count_compar && j != buckets[k])
        {
            k = int(tempValue * (a[j] - mn));
            swap(a[j], a[buckets[k] - 1]);
            --buckets[k];
            ++i;
        }
    }

    for (int i = 1, j; ++count_compar && i < n; ++i) //insertion sort
    {
        j = i;
        while (++count_compar && a[j] < a[j - 1])
        {
            swap(a[j], a[j - 1]);
            --j;
        }
    }
}

//-------------------------------------------------

template <class T>
void HoanVi(T& a, T& b)
{
    T x = a;
    a = b;
    b = x;
}

//-------------------------------------------------

// Hàm phát sinh mảng dữ liệu ngẫu nhiên
void GenerateRandomData(int a[], int n)
{
    srand((unsigned int)time(NULL));

    for (int i = 0; i < n; i++)
    {
        a[i] = rand() % n;
    }
}

void GenerateSortedData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
}

void GenerateReverseData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = n - 1 - i;
    }
}

void GenerateNearlySortedData(int a[], int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
    }
    srand((unsigned int)time(NULL));
    for (int i = 0; i < 10; i++)
    {
        int r1 = rand() % n;
        int r2 = rand() % n;
        HoanVi(a[r1], a[r2]);
    }
}

void GenerateData(int a[], int n, int dataType)
{
    switch (dataType)
    {
    case 0:	// ngẫu nhiên
        GenerateRandomData(a, n);
        break;
    case 1:	// có thứ tự
        GenerateSortedData(a, n);
        break;
    case 2:	// có thứ tự ngược
        GenerateReverseData(a, n);
        break;
    case 3:	// gần như có thứ tự
        GenerateNearlySortedData(a, n);
        break;
    default:
        printf("Error: unknown data type!\n");
    }
}

void runTime(string name, int a[], int n)
{
    cout << fixed << setprecision(5);
    if (name == "selection-sort")
    {
        cout << SelectionSort(a, n) << "s";
        return;
    }
    else if (name == "insertion-sort")
    {
        cout << InsertionSort(a, n) << "s";;
        return;
    }
    else if (name == "bubble-sort")
    {
        cout << BubbleSort(a, n) << "s";;
        return;
    }
    else if (name == "heap-sort")
    {
        cout << HeapSort(a, n) << "s";;
        return;
    }
    else if (name == "merge-sort")
    {
        cout << MergeSort(a, 0, n - 1) << "s";;
        return;
    }
    else if (name == "quick-sort")
    {
        cout << QuickSort(a, 0, n - 1) << "s";;
        return;
    }
    else if (name == "radix-sort")
    {
        cout << RadixSort(a, 0, n - 1, 30) << "s";;
        return;
    }
    else if (name == "shaker-sort")
    {
        cout << ShakerSort(a, n) << "s";;
        return;
    }
    else if (name == "shell-sort")
    {
        cout << ShellSort(a, n) << "s";;
        return;
    }
    else if (name == "counting-sort")
    {
        cout << CountingSort(a, n) << "s";;
        return;
    }
    else
        cout << FlashSort(a, n) << "s";;
    return;
}

void comparsion(string name, int a[], int n)
{
    unsigned long long int  count_compar = 0;
    if (name == "selection-sort")
    {
        SelectionSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "insertion-sort")
    {
        InsertionSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "bubble-sort")
    {
        BubbleSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "heap-sort")
    {
        HeapSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "merge-sort")
    {
        MergeSort_CP(a, 0, n - 1, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "quick-sort")
    {
        QuickSort_CP(a, 0, n - 1, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "radix-sort")
    {
        RadixSort_CP(a, 0, n - 1, 30, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "shaker-sort")
    {
        ShakerSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "shell-sort")
    {
        ShellSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else if (name == "counting-sort")
    {
        CountingSort_CP(a, n, count_compar);
        cout << count_compar;
        return;
    }
    else
    {
        FlashSort_CP(a, n, count_compar);
        cout << count_compar;
    }
    return;
}

int main(int argc, char* argv[])
{
    string mode = argv[1];
    if (mode == "-a")
    {
        cout << "ALGORITHM MODE" << '\n';
        cout << "Algorithm: " << argv[2] << '\n';
        //commandline2
        if (argc > 5)
        {
            string sort = argv[2];
            string size = argv[3];
            string order = argv[4];
            string parameter = argv[5];   
            cout << "Input size: " << argv[3] << '\n';
            cout << "Input order: " << argv[4] << '\n';
            int n = atoi(argv[3]);
            cout << "------------------------------" << endl;
            static int* a = new int[n];
            int dataType = -1;
            if (order == "-rand") {
                dataType = 0;
            }
            else if (order == "-sorted") {
                dataType = 1;
            }
            else if (order == "-reverse") {
                dataType = 2;
            }
            else if (order == "-nsorted") {
                dataType = 3;
            }
            GenerateData(a, n, dataType);
            if (parameter == "-time") {
                cout << "Running time (if required) : ";
                runTime(sort, a, n);
            }
            else if (parameter == "-comp") {
                cout << "Comparisions (if required) : ";
                comparsion(sort, a, n);
            }
            else if (parameter == "-both") {
                cout << "Running time (if required) : ";
                runTime(sort, a, n);
                cout << '\n' << "Comparisions (if required) : ";
                comparsion(sort, a, n);
            }
        }
        else
        {
            string sort = argv[2];
            int n = atoi(argv[3]);
            char* inputfile = argv[3];
            string param = argv[4];
            //commandline3
            if (n != 0)
            {
                cout << "Input size: " << argv[3] << endl;
                cout << endl;
                ofstream f1;
                f1.open("input_1.txt", ios::out | ios::trunc);
                f1 << n << endl;
                cout << "Input order: Randomize" << endl;
                cout << "------------------------------" << endl;
                if (param == "-time") {
                    GenerateRandomData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f1 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-comp") {
                    GenerateRandomData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f1 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-both") {
                    GenerateRandomData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f1 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                f1.close();
                ofstream f2;
                f2.open("input_2.txt", ios::out | ios::trunc);
                f2 << n << endl;
                cout << "Input order: Nearly Sorted" << endl;
                cout << "------------------------------" << endl;
                if (param == "-time") {
                    GenerateNearlySortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f2 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-comp") {
                    GenerateNearlySortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f2 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-both") {
                    GenerateNearlySortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f2 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                f2.close();
                ofstream f3;
                f3.open("input_3.txt", ios::out | ios::trunc);
                f3 << n << endl;
                cout << "Input order: Sorted" << endl;
                cout << "------------------------------" << endl;
                if (param == "-time") {
                    GenerateSortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f3 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-comp") {
                    GenerateSortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f3 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-both") {
                    GenerateSortedData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f3 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                f3.close();
                ofstream f4;
                f4.open("input_4.txt", ios::out | ios::trunc);
                f4 << n << endl;
                cout << "Input order: Reversed" << endl;
                cout << "------------------------------" << endl;
                if (param == "-time") {
                    GenerateReverseData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f4 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-comp") {
                    GenerateReverseData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f4 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                else if (param == "-both") {
                    GenerateReverseData(a, n);
                    for (int i = 0; i < n; i++)
                    {
                        f4 << a[i] << " ";
                        temp[i] = a[i];
                    }
                    cout << "Running time: ";
                    runTime(sort, temp, n);
                    cout << endl;
                    cout << "Comparisions: ";
                    comparsion(sort, temp, n);
                    cout << endl;
                }
                f4.close();
            }
            //commandline1
            else
            {
            int n;
            cout << "Input file: " << inputfile << '\n';
            ifstream intfile;
            intfile.open(inputfile, ios::in);
            intfile >> n;
            cout << "Input size: " << n << '\n';
            cout << "------------------------------" << endl;
            for (int i = 0; i < n; ++i)
            {
                intfile >> a[i];
            }
            static int* a = new int[n];
            if (param == "-time") {
                cout << "Running time (if required) : ";
                runTime(sort, a, n);
            }
            else if (param == "-comp") {
                cout << "Comparisions (if required) : ";
                comparsion(sort, a, n);
            }
            else if (param == "-both") {
                cout << "Running time (if required) : ";
                runTime(sort, a, n);
                cout << '\n' << "Comparisions (if required) : ";
                comparsion(sort, a, n);
            }
            intfile.close();
            }
        }
    }

    //commandline4, 5
    if (mode == "-c")
    {
        cout << "COMPARE MODE" << '\n';
        cout << "Algorithm: " << argv[2] << " | " << argv[3] << '\n';
        string sort1 = argv[2], sort2 = argv[3];
        if (argc == 5)
        {
            int n;
            cout << "Input file: " << argv[4] << '\n';
            char* inputfile = argv[4];
            ifstream intfile;
            intfile.open(inputfile, ios::in);
            intfile >> n;
            for (int i = 0; i < n; ++i)
            {
                intfile >> a[i];
                temp[i] = a[i];
            }

            cout << "Input size: " << n;
            cout << "\n------------------------------\n";

            cout << "Running time: ";
            runTime(sort1, temp, n);
            cout << " | ";
            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            runTime(sort2, temp, n);
            cout << '\n';

            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            cout << "Comparsion: ";
            comparsion(sort1, temp, n);
            cout << " | ";
            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            comparsion(sort2, temp, n);
            intfile.close();
        }
        else
        {
            int n;
            ofstream outfile;
            outfile.open("input.txt", ios::out | ios::trunc);
            cout << "Input file: " << argv[4] << '\n';
            n = atoi(argv[4]);
            outfile << n << '\n';
            cout << "Input oder: " << argv[5];
            cout << "\n------------------------------\n";
            string dataOder = argv[5];
            int dataType = -1;
            if (dataOder == "-rand")
            {
                dataType = 0;
            }
            else if (dataOder == "-nsorted")
            {
                dataType = 3;
            }
            else if (dataOder == "-sorted")
            {
                dataType = 1;
            }
            else
            {
                dataType = 2;
            }

            GenerateData(a, n, dataType);

            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
                outfile << a[i] << ' ';
            }

            cout << "Running time: ";
            runTime(sort1, temp, n);
            cout << " | ";
            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            runTime(sort2, temp, n);
            cout << '\n';

            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            cout << "Comparsion: ";
            comparsion(sort1, temp, n);
            cout << " | ";
            for (int i = 0; i < n; i++)
            {
                temp[i] = a[i];
            }
            comparsion(sort2, temp, n);
            outfile.close();
        }
    }
    return 0;
}