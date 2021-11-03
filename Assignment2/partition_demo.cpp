#include <time.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

// partition a[len] into a[0, 1, ... j-1] <= pivot < a[j, ... len-1], return j
int partition(int a[], int len, int pivot) {
    int low = 0, hi = len-1;
    
    while(1) {
        while (a[low] <= pivot) {
            low++;
            if (low == len) break;
        }
        while (a[hi] > pivot) {
            hi--;
            if (hi < 0) break;
        }

        if (low >= hi) break;
        
        int tmp = a[low];
        a[low] = a[hi];
        a[hi] = tmp;
    }
    return low;
}


int main() {
    int len = 0, pivot = 0;
    int * p;
    int res = partition(p, len, pivot);
    
    cout << res << endl;    
    for (int i = 0; i < len; i++) {
        cout << p[i] << " ";
    }
    cout << endl;
    cout << res << endl;
    
}

