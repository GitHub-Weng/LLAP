//
//  LogUtil.h
//  llap
//
//  Created by wengdada on 03/02/2018.
//  Copyright © 2018 Shenzhen University. All rights reserved.
//

#ifndef LogUtil_h
#define LogUtil_h
#define TURNOFFLOG 0
#include <cstdio>
#define TRACE_FUN(...) \
printf("wengdada call function <%s>: \n",__FUNCTION__, ##__VA_ARGS__)
inline void log1VFloat32(float* vector,int num)
{
    TRACE_FUN();
    if(TURNOFFLOG)
    {
        return;
    }
    
    float  temp = 0;
    int i = 0;
    for(i= 0;i<num;i++)
    {
        temp = vector[i];
        printf("%f ",temp);
    }
}

inline void log2VFloat32(float* vector,int num1, int num2)
{
    TRACE_FUN();
    if(TURNOFFLOG){
        return;
    }
    float  temp = 0;
    int i = 0;
    for(i= 0;i<num1;i++)
    {
        int j = 0;
        for(j=0;j<num2;j++)
        {
            temp = vector[i];
            printf("%f ",temp);
        }
    }
}

inline void log3VFloat32(float* vector,int num1, int num2,int num3)
{
    TRACE_FUN();
    if(TURNOFFLOG)
    {
        return;
    }
    float  temp = 0;
    int i = 0;
    for(i= 0;i<num1;i++)
    {
        int j = 0;
        for(j=0;j<num2;j++)
        {
            int z = 0;
            for(z=0;z < num3;z++)
            {
                temp = vector[i];
                printf("%f ",temp);
            }
        }
    }
}


inline void log1VInt16(short* vector,int num)
{
    TRACE_FUN();
    if(TURNOFFLOG)
    {
        return;
    }
    short  temp = 0;
    int i = 0;
    for(i= 0;i<num;i++)
    {
        temp = vector[i];
        printf("%i ",temp);
    }
}



#endif /* LogUtil_h */

