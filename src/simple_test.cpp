#include <iostream>
#include <string>
#include <vector>
using namespace std;

int main()
{
    vector<int> v{1, 2, 3, 4, 5, 6, 7, 8};
    for (auto &i : v)
        i *= i;
    for (auto i : v)
        cout << i << " ";
    cout << endl;
    cout << v.empty() << endl
         << v.size() << endl;
}