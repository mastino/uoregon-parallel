// While not particulary secure... it turns out this strategy also isn't too 
// bad (as long as the keys are kept secret and have different lengths and 
// aren't reused together)... or at least that is what a member of the security
// group that works in crypto told me.

#include "key.h"
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include <tbb/tbb.h>
#include <pthread.h>

using namespace tbb;
using namespace std;

// utility function: given a list of keys, a list of files to pull them from, 
// and the number of keys -> pull the keys out of the files, allocating memory 
// as needed
void getKeys(xorKey* keyList,char** fileList, int numKeys) {
  int keyLoop=0;
  for(keyLoop=0;keyLoop<numKeys;keyLoop++) {
     readKey(&(keyList[keyLoop]), fileList[keyLoop]);
  }
}

// A function to perform the basic reduction operation.
// here two arguments are reduced to one result of the same type
char whatever(char a, char b) {
		char c = a ^ b ; 
		return c;
}

// A pointer to the reducing operator to pass over the collection 
char (*functionPtr)(char,char) = whatever;

// Given text, a list of keys, the length of the text, 
// and the number of keys, encodes the text with the supplied function
void encode(
	char (*functionPtr)(char,char), //the supplied function
	char* plainText, 
	char* cypherText, 
	xorKey* keyList, 
	int ptextlen, 
	int numKeys
) {

  
  task_scheduler_init (10);
  parallel_for( blocked_range<int>(0, ptextlen),
    [&]( const blocked_range<int> &r ) { 

        for(int charLoop = r.begin(), el_end = r.end(); charLoop <= el_end; charLoop++){ 

          int keyLoop=0;
          char cipherChar=plainText[charLoop]; 
          for(keyLoop = 0; keyLoop < numKeys; keyLoop++) {
             cipherChar = (*functionPtr) (cipherChar, getBit( &(keyList[keyLoop]), charLoop));
          }

          cypherText[charLoop]=cipherChar;
        } 

    });


}

void decode(char (*functionPtr)(char,char),
	char* cypherText, 
	char* plainText, 
	xorKey* keyList, 
	int ptextlen, 
	int numKeys
) {
  encode(functionPtr,cypherText, plainText, keyList, ptextlen, numKeys); 
  //isn't symmetric key cryptography awesome? 
}

int main(int argc, char* argv[]) {
  if(argc<=2){
      printf("Usage: %s <fileToEncrypt> <key1> <key2> ... <key_n>\n",argv[0]);
      return 1;
  }

  // read in the keys
  int numKeys=argc-2;
  xorKey* keyList=(xorKey*)malloc(sizeof(xorKey)*numKeys); // allocate key list
  getKeys(keyList,&(argv[2]),numKeys);
  
  // read in the data to encrypt/decrypt
  off_t textLength=fsize(argv[1]); //length of our text
  FILE* rawFile=(FILE*)fopen(argv[1],"rb"); //The intel in plaintext
  char* rawData = (char*)malloc(sizeof(char)*textLength);
  fread(rawData,textLength,1,rawFile);
  fclose(rawFile);

  // Encrypt
  char* cypherText = (char*)malloc(sizeof(char)*textLength);
  encode(functionPtr,rawData,cypherText,keyList,textLength,numKeys);

  // Decrypt
  char* plainText = (char*)malloc(sizeof(char)*textLength);
  decode(functionPtr,cypherText,plainText,keyList,textLength,numKeys);

  // Check
  int err = 0;
  int i;
  for(i=0;i<textLength;i++) {
    if(rawData[i]!=plainText[i]) {
      fprintf(stderr, "Encryption/Decryption is non-deterministic\n");
      i=textLength;
      err = 1;
    }
  }
 
  if (err == 0) 
    fwrite (cypherText, sizeof(char), textLength, stdout);
  
  return err;

}
