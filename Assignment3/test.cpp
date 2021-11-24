#include <limits.h>
#include <iostream>
#include <limits>
using namespace std;

int main() {
    cout << INT_MAX << endl;
    cout << INT_MAX + 1 << endl;

    cout << numeric_limits<double>::infinity() << endl;
    cout << ((numeric_limits<double>::infinity()  + 100) > 1000) << endl;

    return 0;
}