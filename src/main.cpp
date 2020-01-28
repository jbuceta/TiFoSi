#include <iostream>
using namespace std;

#include "auxiliar.h"
#include "main.h"

int main()
{
  poblacio p;
  
  p.crea_poblacio();
  
  p.bucle_dinamica_principal();
  
  p.destrueix_poblacio();
  
  return 0;
}
