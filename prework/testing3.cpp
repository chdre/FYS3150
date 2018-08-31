#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

int main(){

  int greetingOption = 2;

  switch(greetingOption){

    case 1:
      cout << "bonjour" << endl;
      break;
    case 2:
      cout << "hola" << endl;
      break;
    default :
    cout << "hello" << endl;
  }

  return 0;
}
