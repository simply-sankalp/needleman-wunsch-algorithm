//  dependencies

#include <stdio.h>
#include <string.h>
#define MAX_LENGTH 500
#define gap_penalty -10
#define extend_penalty -0.5

// function declarations

int match_mismatch(char a, char b);
int max(int a, int b, int c);
int len(char str[]);
int similarity(char align_A[], char align_B[]);
float score(char align_A[], char align_B[]);

// Needleman Wunsch Algorithm

void needleman_wunsch(char seq_A[], char seq_B[])
{

    //  forming the score matrix

    int m = len(seq_A) + 1,
        n = len(seq_B) + 1;

    int matrix[m][n];

    //  initializing the matrix

    matrix[0][0] = 0;

    for (int i = 1; i < m; i++)
    {
        matrix[i][0] = gap_penalty * i;
    }

    for (int j = 1; j < n; j++)
    {
        matrix[0][j] = gap_penalty * j;
    }

    //  scoring

    for (int i = 1; i < m; i++)
    {
        for (int j = 1; j < n; j++)
        {
            int case_1 = matrix[i-1][j-1] + match_mismatch(seq_A[i-1], seq_B[j-1]),
                case_2 = matrix[i-1][j] + gap_penalty,
                case_3 = matrix[i][j-1] + gap_penalty;
            
            matrix[i][j] = max(case_1, case_2, case_3);
        }
    }

    // for (int i = 0; i < m; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         printf("%d\t", matrix[i][j]);
    //     }
    //     printf("\n");
    // }

    //  finding the best route along the matrix

    int x = m-1,
        y = n-1;
    char ch;
    
    char movements[MAX_LENGTH * 2];
    movements[0] = '\0';

    //printf("\n( %d , %d )\n", x, y);
    
    while (x > 0 || y > 0)
    {
        if (x == 0)
        {
            ch = 'T';
        }
        else if (max(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) == matrix[x-1][y-1])
        {
            ch = 'D';
            x--;
            y--;

            //printf("( %d , %d )\n", x, y);
        }
        else if (max(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) == matrix[x][y-1])
        {
            ch = 'L';
            y--;

            //printf("( %d , %d )\n", x, y);
        }
        else
        {
            ch = 'T';
            x--;

            //printf("( %d , %d )\n", x, y);
        }
        strncat(movements, &ch, 1);
    }

    int size = len(movements);

    // printf("\n\n");

    // for (int i = 0; i < size; i++)
    // {
    //     printf("%c", movements[i]);
    // }

    // printf("\n\n");


    //  storing the aligned sequences

    char align_A[size+1];
    char align_B[size+1];
    align_A[0] = '\0';
    align_B[0] = '\0';

    int ind_A = m-2,
        ind_B = n-2;
    
    char ch_A, ch_B;

    for (int i = 0; i < size; i++)
    {
        if (movements[i] == 'D')
        {
            ch_A = seq_A[ind_A];
            ch_B = seq_B[ind_B];
            ind_A--;
            ind_B--;
        }
        else if (movements[i] == 'L')
        {
            ch_A = '-';
            ch_B = seq_B[ind_B];
            ind_B--;
        }
        else
        {
            ch_A = seq_A[ind_A];
            ch_B = '-';
            ind_A--;
        }

        strncat(align_A, &ch_A, 1);
        strncat(align_B, &ch_B, 1);
    }

    //  computing gap

    int gap = 0;

    for (int i = 0; i < size; i++)
    {
        if (align_A[i] == '-' || align_B[i] == '-')
        {
            gap++;
        }
    }

    //  computing identity
    
    int identity = 0;
    
    for (int i = 0; i < size; i++)
    {
        if (align_A[i] == align_B[i])
        {
            identity++;
        }
    }

    float identity_percentage = ((float)identity/size)*100,
          similarity_percentage = ((float)similarity(align_A, align_B)/size*100),
          gap_percentage = ((float)gap/size)*100;
    

    //  output

    printf("\n#\tGap Penalty :\t%d", gap_penalty);
    printf("\n#\tExt Penalty :\t%.1f", extend_penalty);
    printf("\n#\tLength :\t%d", len(movements));
    printf("\n#\tIdentity :\t%d/%d  (%.2f%%)", identity, size, identity_percentage);
    printf("\n#\tSimilarity :\t%d/%d  (%.2f%%)", similarity(align_A, align_B), size, similarity_percentage);
    printf("\n#\tGaps :\t\t%d/%d  (%.2f%%)", gap, size, gap_percentage);
    printf("\n#\tScore :\t\t%.1f", score(align_A, align_B));
    printf("\n#");
    printf("\n#\t");

    //printing aligned sequence A

    for (int i = 1; i <= len(align_A); i++)
    {
        printf("%c", align_A[len(align_A) - i]);
    }


    //printing the connectors

    printf("\n#\t");

    for (int i = size - 1; i >= 0; i--)
    {
        if (align_A[i] == align_B[i])
        {
            printf("|");
        }
        else if (align_A[i] == 'A' && align_B[i] == 'G')
        {
            printf(".");
        }
        else if (align_A[i] == 'G' && align_B[i] == 'A')
        {
            printf(".");
        }
        else if (align_A[i] == 'C' && align_B[i] == 'T')
        {
            printf(".");
        }
        else if (align_A[i] == 'T' && align_B[i] == 'C')
        {
            printf(".");
        }
        else
        {
            printf(" ");
        }
    }

    //printing aligned sequence B

    printf("\n#\t");

    for (int i = 1; i <= len(align_B); i++)
    {
        printf("%c", align_B[len(align_B) - i]);
    }

    printf("\n");
}

// driver code

int main()
{
    //  taking input
    
    char filename1[100];
    char filename2[100];

    printf("\nPlease enter address of first file: ");
    scanf("%s", filename1);
    printf("\nPlease enter address of second file: ");
    scanf("%s", filename2);

    char seq_A[MAX_LENGTH];
    char seq_B[MAX_LENGTH];

    FILE *file1 = fopen(filename1, "r");
    if (file1 == NULL) {
        printf("Error opening file %s\n", filename1);
        return 1;
    }

    if (fgets(seq_A, MAX_LENGTH, file1) == NULL) {
        printf("Error reading string from file %s\n", filename1);
        fclose(file1);
        return 1;
    }

    fclose(file1);

    FILE *file2 = fopen(filename2, "r");
    if (file2 == NULL) {
        printf("Error opening file %s\n", filename2);
        return 1;
    }

    if (fgets(seq_B, MAX_LENGTH, file2) == NULL) {
        printf("Error reading string from file %s\n", filename2);
        fclose(file2);
        return 1;
    }

    fclose(file2);
    
    //  output

    printf("\n\n");
    for (int i = 0; i < 31; i++)
    {
        printf("#");
    }

    printf("\n#");
    printf("\n#\tAligned Sequences : 2");
    printf("\n#\t1 : %s", filename1);
    printf("\n#\t2 : %s", filename2);

    printf("\n#");

    needleman_wunsch(seq_A, seq_B);

    printf("#\n");
    for (int i = 0; i < 31; i++)
    {
        printf("#");
    }

    printf("\n\n");

    return 0;
}

//function definitions

int match_mismatch(char a, char b)
{
    //adheres to EDNAFULL matrix

    if (a == b) return 5;
    else return -4;
}

int max(int a, int b, int c)
{
    if (a >= b && a >= c)
    {
        return a;
    }
    else if (b >= c && b >= a)
    {
        return b;
    }
    else return c;
}

int len(char str[])
{
    int i = 0;
    while(str[i] != '\0')
    {
        i++;
    }
    return i;
}

int similarity(char align_A[], char align_B[])
{
    int val = 0;

    for (int i = 0; i < len(align_A); i++)
    {
        if (align_A[i] == 'A' || align_A[i] == 'G')
        {
            if (align_B[i] == 'A' || align_B[i] == 'G')
            {
                val++;
            }
        }
        else if (align_A[i] == 'C' || align_A[i] == 'T')
        {
            if (align_B[i] == 'C' || align_B[i] == 'T')
            {
                val++;
            }
        }
        else if (align_A[i] == '-' && align_B[i] == '-')
        {
            val++;
        }
    }

    return val;
}

float score(char align_A[], char align_B[])
{
    float val = 0;

    for (int i = 0; i < len(align_A); i++)
    {
        if (align_A[i] == '-')
        {
            if (i == 0 || align_A[i-1] != '-')
            {
                val += gap_penalty;
            }
            else if (i != 0 && align_A[i-1] == '-')
            {
                val += extend_penalty;
            }
            else
            {
                val += match_mismatch(align_A[i], align_B[i]);
            }
        }
    }

    return val;
}