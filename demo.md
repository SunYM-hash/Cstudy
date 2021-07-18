### C++笔记

***

#### 一、C++ primer 阅读

1、内置变量在函数体外定义时初始化为0（包括main函数内），在外不自动初始化，会生成随机值。

2、const指向的变量为文件局部变量，想被其他文件应用需加入extern声明。

```c++
extern const int a=0;    // const变量定义时必须初始化
```

3、引用必须关联到其他变量，引用不能改变对象，const变量应用必须const引用。const引用可以绑定右值或者相关类型但不同类型的变量。

``` c++
int & a;         // error
int & a = 10;  //error
double b;
const int & a = 10, &c = b; // 合法
```

4、指针与数组区分

``` c++
int *p[n]    //   指针数组
int (*p)[n]  // 指向一维数组的指针
```

5、位运算符

| 符号  |    操作    |
| :---: | :--------: |
|   &   |     与     |
|  \|   |     或     |
|   ^   |    异或    |
|  ～   |    取反    |
| <<,>> | 左移和右移 |

6、指向const对象的指针和const指针（要理解）。

```c++
const int *p ;    //指向const对象的指针
int *const p = &b ;   // const指针,必须初始化
```

指向const对象的指针，不能用其来改变其指向对象的值，所以可以将形参定义为指向const对象的指针，避免修改实参。（这里指向const对象的指针指向的不一定是const对象，实际上编译器并不知道你指向的是不是const对象，所以默认是const对象）。同时要注意，因为const指针指向的不一定是const对象，所以如果对象不是const对象，那对象虽然不能被指针修改，但可以被其他方式修改。

7、C风格字符串，一定要有"\0"。在编写中尽量避免使用c风格字符串，改用string既提高效率，也更安全。strlen不包括空字符，使用cstring库时，需计算大小。

8、new分配类型时，如果为类对象，使用该类对象的默认构造函数。如果为内置类型，则无初始化。

``` c++
int* p = new int[10];   //不初始化
int* p = new int[10]();//初始化为0，（）要求编译器进行初始化。
int* p = new string[10]; //已初始化

delete [] p; // 对于动态数组，如果漏掉方括号，编译器会认为指向单个对象，导致内存泄漏。  
```

9、如果引用的过程中避免复制实参是唯一目的，应将形参定义为const引用

``` c++
bool isShoter(const string &s1,const string & s2)
```

10、const引用，只读不改。用引用做形参时，最好避免普通引用做实参，而是改用const引用。

``` c++
int icur(int &val){
 	retrun icur++;
}
int main(){
    const int v1 = 10;
	int v3 = incr(v1); //error
    int v4 = incr(0);  // error
    
}
```



***



#### 二、LeetCode刷题心得

1、C++位运算符，<<,>>

``` c++
int mid = left + (right-left)/2;   //二分查找
int mid = left +((right-left)>>1);   //位运算，速度更快
if (a %2 == 1)   if (a&1)   // 与运算，判断奇偶
```

2、数组增删元素复杂度均为O(n)

3、if，else if有先后判断顺序，记得先判断边界条件

4、map中都是键值对，只能通过键找值，不能靠值来找到键（得自己实现）。

5、[最小回文子串](https://leetcode-cn.com/problems/longest-palindromic-substring/) 题可得，动态规划题要注意根据初始化的状态来决定动态规划的遍历顺序。(能一起判断尽量一起判断，不然可能产生不好发现的bug)

|               | dp[i] [j] |
| :-----------: | :-------: |
| dp[i+1] [j-1] |           |

6、typeid(i).name()查看变量类型，必须包含头文件<typeinfo> 

***



#### 三、Github资料

1、常用命令

``` bash
git init   #初始化仓库
git add 1.txt   #提交到缓冲区
git commit -m "备注"　　　#将缓冲区的东西提交到归档去
git remote add origin http..    #添加远程仓库
git push origin master     #提交到远程仓库
```











