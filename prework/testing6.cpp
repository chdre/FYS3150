#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

int main(){

  int myFavNums[5];

  int badNumber[5] = {4, 13, 14, 24, 34};

  cout << "Bad number 1: " << badNumber[0] << endl;

  char myName[5][5] = {{'C', 'h', 'r', 'i', 's'}, {'D', 'r', 'e', 'i'}};

  cout << "2nd letter 2nd arr " << myName[1][1];

  myName[0][2] = 'e';

  cout << "New value " << myName[0][2] << endl;

  for(int i = 1; i <= 10; i++){

    cout << i << endl;
  }

  for(int j = 0; j < 2; j++){

    for(int k = 0; k < 5; k++){

      cout << myName[j][k];
    }
    cout << endl;
  }

  return 0;
}
