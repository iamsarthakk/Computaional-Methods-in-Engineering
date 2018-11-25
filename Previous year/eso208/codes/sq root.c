#include <stdio.h>

int main(void) {
	float n;
	scanf("%f",&n);
	int i;
	float j=1,k;
	for(i=0;i<100;i++)
	{
	    k = (j + n/j)*0.5;
	    j = k;
	}
	printf("%f",j);
	return 0;
}

