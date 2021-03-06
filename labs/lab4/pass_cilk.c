#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#include <cilk/cilk.h>
#include <pthread.h>
#include <unistd.h>

#include <openssl/md5.h>

const char* chars="0123456789";

// tests if a hash matches a candidate password
int test(const char* passhash, const char* passcandidate) {
    unsigned char digest[MD5_DIGEST_LENGTH];
    
    MD5((unsigned char*)passcandidate, strlen(passcandidate), digest);
    
    char mdString[34];
    mdString[33]='\0';
    for(int i=0; i<16; i++) {
        sprintf(&mdString[i*2], "%02x", (unsigned char)digest[i]);
    }
    
    return strncmp(passhash, mdString, strlen(passhash));
}

// maps a PIN to a string
void genpass(long passnum, char* passbuff) {
    passbuff[8]='\0';
    int charidx;
    int symcount=strlen(chars);
    for(int i=7; i>=0; i--) {
        charidx=passnum%symcount;
        passnum=passnum/symcount;
        passbuff[i]=chars[charidx];
    }
}

/* Return the password if matched successfully, NULL otherwise */
char * password(char * passmatch, char * digest, long element)
{
   genpass(element, passmatch);
   if (test(digest, passmatch) == 0) {
      return passmatch;
   }

   return NULL;
}

// generates potential passwords for a given hash value and returns when one is found
// params fp = a function pointer to generate the password
//        passmatch = memory buffer used in computing password
//        digest = the hash to be located
//        max_val = the bggest value to try?
char * map_reduce( char* (fp)(char*,char*,long), char * passmatch, char * digest, long max_val) {

  long element;
  int i;

  char done = 0;
  pthread_mutex_t d;


  pthread_mutex_init(&d, NULL);
  cilk_for (element = 0; element <= max_val; element++) {
    if(done == 0) {
      char buffer[9];
      char* match;
      match = fp(buffer, digest, element);
      if ((done == 0) && (match != NULL)) {
        pthread_mutex_lock(&d);
        for(i = 0; i < 9; i++) passmatch[i] = match[i];
        done = 1;
        pthread_mutex_unlock(&d);
      }
    }
  }


  if(done != 0)
    return passmatch;
  else
    return NULL;
}

int main(int argc, char** argv)
{
    char passmatch[9], * match;

    if (argc != 2) {
        printf("Usage: %s <password hash>\n",argv[0]);
        return 1;
    }

    match = map_reduce(password, passmatch, argv[1], 99999999);

    if (match == NULL) {
       printf("\nERROR password not found for: %s\n\n", argv[0]);
    } else {
       printf("\nSUCCESS found: %s\n\n", match);
    }

    return 0;
}
