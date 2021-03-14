#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <queue>
#include <algorithm>
//pentru timp in milisecunde
//#include <chrono>
using namespace std;

ifstream f("file.txt");

clock_t start, endt, execution_time;

void print_vector(vector <int> &v)
{
    for(int i = 0; i < v.size(); ++i)
        cout<<v[i]<<" ";

    cout<<endl;
}



bool test_sort(vector <int> copie, vector <int> &v)
{
    vector <int> frec(copie.size(), 0);

    if(copie.size() != v.size())
        return false;

    //verificam daca valorile sunt ascendente si retinem in vectorul de frecventa valorile din copie
     for(int i = 0; i < copie.size(); i++)
    {
        if(i != copie.size() - 1 && copie[i] > copie[i+1])
        {
            return false;}

        frec[copie[i]]++;
    }

    //eliminam valorile din vector
    for(int i = 0; i < v.size(); i++)
        frec[v[i]]--;

    for(int i = 0; i < frec.size(); i++)
        if(frec[i] != 0)
            return false;

    return true;
}


void bubble_sort(vector <int> &v, int n)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n - i - 1; j++)
        {
            if (v[j] > v[j + 1])
            {
                int aux = v[j];
                v[j] = v[j+1];
                v[j+1] = aux;
            }
        }
    }
}

void count_sort(vector <int> &v, int n, int maxim){
    vector <int> fr(maxim + 1, 0);
    vector <int> sorted;
    for(int i = 0; i < n; i++)
    {
        fr[v[i]]++;
    }
    for(int i = 0; i <= maxim; i++)
    {
        while(fr[i] != 0)
        {
            sorted.push_back(i);
            fr[i]--;
        }
    }
    v = sorted;
}

void merge(vector <int> &v, int left, int middle, int right)
{
    int i = left, j = middle + 1, k = 1;

    vector <int> aux;

    while(j <= right && i <= middle){
        if(v[j] > v[i]){
            aux.push_back(v[i++]);
        }
        else{
            aux.push_back(v[j++]);
        }
    }

    while(i <= middle){
        aux.push_back(v[i++]);
    }
    while(j <= right){
        aux.push_back(v[j++]);
    }

    for(i = left; i <= right; i++){
        v[i] = aux[i-left];
    }

}

void merge_sort(vector <int> &v, int left, int right)
{
    if(right > left){
        int middle = (left + right) / 2;
        merge_sort(v, left, middle);
        merge_sort(v, middle + 1, right);

        merge(v, left, middle, right);
    }
}

void swap_values(vector <int> &v, int left, int right)
{
    int aux = v[left];
    v[left] = v[right];
    v[right] = aux;
}

void swap(int &a, int &b)
{
    int aux = a;
    a = b;
    b = aux;
}

int partition_by_pivot( vector <int> &v, int left, int right, int pivot)
{
    while(left <= right)
    {
        while(v[left] < pivot) left++;
        while(v[right] > pivot) right--;

        if(left <= right)
        {
            int aux = v[left];
            v[left] = v[right];
            v[right] = aux;
            left++;
            right--;
        }
    }
    return left;
}

void quicksort_med_3(vector <int> &v, int left, int right)
{
    if(left < right)
    {
        int middle = (left + right) / 2, pivot;

        if(v[left] > v[right])
            swap_values(v, left, right);

        if(v[left] > v[middle])
            swap_values(v, middle, left);

        if(v[middle] > v[right])
            swap_values(v, right, middle);

        pivot = v[middle];
        int index = partition_by_pivot(v, left, right, pivot);

        quicksort_med_3(v, left, index-1);
        quicksort_med_3(v, index, right);
    }
}


void quick_sort_random(vector <int> &v, int left, int right)
{
    int l = left, r = right; // i si j
    int pivot_rand = v[rand() % (right - left + 1) + left];

    while( l < r)
    {
        while(v[l] < pivot_rand) l++;
        while(v[r] > pivot_rand) r--;

        if( l <= r)
        {
            swap(v[l++],v[r--]);
        }
    }

    if(r > left)
        quick_sort_random(v, left, r);
    if(l < right)
        quick_sort_random(v, l, right);
}

int maximum_value(vector <int> &v, int N)
{
    int max_value = v[0];
    for (int i = 1; i < N; i++)
        if (v[i] > max_value)
            max_value = v[i];
    return max_value;
}

void count_sort_radix(vector <int> &v, int N, int exp)
{
    vector <int> ans(N,0);
    vector <int> count(N,0);

    for(int i = 0; i < N; ++i)
        count[(v[i] / exp) %10] ++ ;

    for(int i = 1; i < 10; ++i)
        count[i] += count[i-1];

    for(int i = N-1; i >= 0; i--)
    {
        ans[count[(v[i] / exp) %10] - 1] = v[i];
        count[ ( v[i] / exp ) % 10]--;
    }

    v = ans;

}

vector <int> radix_sort(vector <int> v, int N, int Max)
{
    int mx = maximum_value(v, N);

    for( int exp = 1; mx/exp > 0; exp *=10)
        count_sort_radix(v, N, exp);

    return v;
}

void radix_sort_base_16(vector <int> &v, int N)
{
    //lista de cozi
    queue <int> q[16];
    int k;

    for( int i = 0; i <= 24; i += 4 )
    {
        for(int j = 0; j < N; ++j)
            q[( v[j] >> i)& 15].push(v[j]);

        k = 0;

        for(int j = 0; j <= 15 && k <  N; ++j)
            while(!q[j].empty())
            {
                v[k++] = q[j].front();
                q[j].pop();
            }
    }

}

void radix_sort_base_4(vector <int> &v, int N)
{
    queue<int> q[256];
    int k;

    for(int i = 0; i <= 24; i += 2)
    {
        for(int j = 0; j < N; ++j)
            q[(v[j] >> i) & 3].push(v[j]);

        k = 0;

        for(int j = 0; j <= 3 && k < N; ++j)
            while(!q[j].empty())
            {
                v[k++] = q[j].front();
                q[j].pop();
            }
    }
}


int main()
{
    int T, N, Max;
    vector <int> v1; //v1 va fi vectorul format din elementele randomizate
    f>>T;

    for(int i = 0; i < T; ++i)
    {
        f>>N>>Max;

        for(int j = 0 ; j < N; j++)
        {
            v1.push_back(rand() % (Max + 1));
        }

        cout<<"N = "<<N<<" Max = "<<Max<<endl;

        vector <int> v_copy = v1;

        //Bubble sort
        cout<<"Pentru Bubble sort: ";
        if(N > 10000){
            cout<<"Vectorul nu poate fi sortat utilizat metoda Bubble sort.";
        }
        else{
            //vector <int> v_copy =  v1;
            //clock_t start_t = clock();
            start = clock();
            bubble_sort(v_copy, N);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
            if(test_sort(v_copy, v1)){
                cout<<"succes";
            }

            /*
            pentru timp in milisecunde

            auto t1 = high_resolution_clock::now();
            bubble_sort(v_copy, N);
            auto t2 = high_resolution_clock::now();
            duration<double, std::milli> ms_double = t2 - t1;
            cout<<"BUBBLE SORT TIME: "<<ms_double.count()<<"ms    ";*/
        }

        //Count sort
        cout<<endl<<"Pentru Count sort: ";
        if(N > 10000){
            cout<<"Vectorul nu poate fi sortat utilizat metoda Count sort.";
        }
        else{
                v_copy.clear();
                v_copy =  v1;
                start = clock();
                count_sort(v_copy, N, Max);
                endt = clock();
                execution_time = double(endt - start)/CLOCKS_PER_SEC;
                cout<<"timp : " <<execution_time<<endl;
                 if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
            }

        //Merge sort
        if(N < 100000000)
        {
            cout<<endl<<"Pentru Merge sort: ";
            v_copy = v1;
            start = clock();
            merge_sort(v_copy,0,N-1);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
           // execution_time = double(endt - start)/CLOCKS_PER_SEC;
           // cout<<"timp : " <<execution_time<<endl;
             if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else{cout<<"Nu se poate efectua Merge sort.";}




        //Quicksort
        //1
        if(N < 100000000)
        {

            cout<<endl<<"Pentru Quick sort cu mediana 3: ";
            v_copy = v1;
            start = clock();
            quicksort_med_3(v_copy, 0, N-1);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
             if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else{
            cout<<"Nu se poate efectua Quicksort cu pivot mediana 3.";
        }


        //2
        if(N < 100000000)
        {
            cout<<endl<<"Pentru Quick sort cu pivot ales random: ";
            v_copy = v1;
            start = clock();
            quick_sort_random(v_copy,0,N-1);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
             if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else
        {
            cout<<"Nu se poate efectua Quicksort cu pivot randomizat.";
        }



        //Radix Sort
        v_copy = v1;
        if(N < 100000000)
        {
            cout<<endl<<"Pentru Radix sort cu baza 10: ";
            start = clock();
            v_copy = radix_sort(v_copy, N, Max);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
            if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else
        {
            cout<<"Nu se poate efectua Radixsort in baza 10.";
        }


        if(N < 100000000)
        {
            cout<<endl<<"Pentru Radix sort cu baza 4: ";
            v_copy = v1;
            start = clock();
            radix_sort_base_4(v_copy,N);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
            if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else{
            cout<<"Nu se poate efectua Radix sort in baza 4.";
        }


        if(N < 100000000)
        {
            cout<<endl<<"Pentru Radix sort cu baza 16: ";
            v_copy = v1;
            start = clock();
            radix_sort_base_16(v_copy,N);
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
            if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
        }
        else
        {
            cout<<"Nu se poate efectua Radix sort in baza 16.";
        }


        if(N < 100000000 && Max < 1000000)
        {
            cout<<endl<<"Pentru Sort: ";
            v_copy = v1;
            start = clock();
            sort(v_copy.begin(), v_copy.end());
            endt = clock();
            execution_time = double(endt - start)/CLOCKS_PER_SEC;
            cout<<"timp : " <<execution_time<<endl;
            if(test_sort(v_copy, v1)){
                    cout<<"succes";
            }
            cout<<test_sort(v_copy, v1);
        }
        else{
            cout<<"Nu se poate efectua sortarea nativa.";
        }

        v1.clear();
    }

    return 0;
}
