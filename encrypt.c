#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//does not work with any characters other than spaces and letters

int main()
{
    //holds the original message
    char sentence[200];

    printf("Enter a sentence\n");
    scanf("%[^\n]", sentence);


    //converts all letters to lowercase
    for(int i=0; i<strlen(sentence); i++){
        if(sentence[i] <= 96){
           sentence[i] += 32;
        }
    }

    //this encrypts the string and the key undoes it
    int cypher = 0;
    int key = 0;

    //encrypt holds the encrypted sentence
    int encrypt[200];
    
    //the number of spaces
    int number_spaces = 0;

    //decrypt holds the encrypted sentence and room for a space after each character
    int decrypt[400];

    //converts the sentence to numbers
    for(int i=0; i<strlen(sentence); i++){
        //increases the cypher by i squared for each letter so that repeated letters are different
        cypher += i*i;
        
        encrypt[i] = sentence[i] + cypher;

	number_spaces++;
    }

    //prints out the encryption with the number of spaces squared plus seven at the beginning
    //this is required for the real decrypt.c code
    printf("%i ", number_spaces*number_spaces+7);
    for(int i=0; i<strlen(sentence); i++){
        printf("%i ", encrypt[i]);
    }


    //the decryption here is a check to make sure it was encrypted properly
    printf("\n");
   
    for(int i=0; i<=number_spaces; i++){
        //undoes the encryption
        key -= i*i;
        decrypt[i] = encrypt[i] + key;

        if(decrypt[i] <= 96){
           decrypt[i] -= 32;
        }

        //undoes the cipher and prints the check
        if(decrypt[i] >= 32 && decrypt[i] <= 122){
            printf("%c", decrypt[i]);
        }
    }

    printf("\n");


    return 0;
}

