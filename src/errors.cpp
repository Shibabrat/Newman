// errors.cpp
//
// error module
//
//   Author: Philip Du Toit 
//   Rev 1.0
//

#include "errors.h"

#include <iostream>

void FatalError(const char error_text[]) {
  cout << error_text << endl;
  cout << "\n Exiting to system ...\a" << endl;
  exit(1);
}

