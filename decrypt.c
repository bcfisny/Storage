#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//does not work with any characters other than spaces and letters

int main()
{

    int number_spaces;

    //the number of spaces need to be counted in the encoded message
    scanf("%i", &number_spaces);

    //the stopping condition of the loop is the number of characters in the message
    int length = sqrt(((number_spaces-7)));

    //encrypt holds the encrypted numbers (the input) and decrypt displays the output
    int *encrypt = (int *)malloc((number_spaces) * sizeof(int));
    int *decrypt = (int *)malloc((number_spaces) * sizeof(int));

    //scans the input and counts the number of characters since there is a space after each one
    for(int i=0; i<length; i++){
        scanf("%i", &encrypt[i]);
    }


    int key = 0;

    printf("\n");
    //decrypts the string
    for(int i=0; i<number_spaces; i++){
        //undoes the encryption
        key -= i*i;
        decrypt[i] = encrypt[i] + key;

        if(decrypt[i] <= 96){
           decrypt[i] -= 32;
        }

        //prints the original message
        if(decrypt[i] >= 32 && decrypt[i] <= 122){
            printf("%c", decrypt[i]);
        }
    }

    printf("\n");

    free(encrypt);
    free(decrypt);

    return 0;
}


