## Note for `pauli_string.py`

### 写在最前面

**此前使用的对每个不同的Z_1，生成一组generator的做法实际上并不能够穷举所有的情况**

这部分代码仍然被保留，写在`Get_Random_Generators`以及`Get_Minimun_Cases`中。

### 主要框架

代码中主要是完成了两个类的构建：
- `Pauli_Operator`：它用于存储每个Pauli算符对应的`sympletic_form`，直积形式（`pauli_string`）以及长度（`length`）。
- `Clifford_Operator`：由于希望便利所有可能的`Clifford_Operator`（虽然最终没有实现），代码中设计了两种不同的生成`Clifford_Operator`的方案：第一种是直接输入Pauli群的生成元经过Clifford变换后的`sympletic_form`。第二种方案是输入`Clifford_Operator`的直积形式，它是一个list，其中每一个元素都是一个tuple，tuple中第一个位置为`"H"`或`"S"`或`"CNOT"`（它们是Clifford group的一组生成元，当然，是不考虑相位的情况），第二个位置是这些gate作用的qubit位置。该class中定义了从第二种方案到第一种方案的转化（但是没有从第一种方案到第二种方案的转化，原因很简单，这个方向我不会hh）。

### 分析讨论

我们将作用在n个量子比特上的Pauli群中原长度为n的Pauli算符的集合记为$\mathcal{P}$。上述Clifford group的生成元中，实际上能够改变集合$\mathcal{P}$的平均长度的只有CNOT，因为H仅仅只是将该位置上的X和Z互换，S仅仅只是将X和Y互换，因此H和S单独作用并不能改变$\mathcal{P}$的平均长度。

#### Result 0

各个Clifford generator的作用：
```
For Hadamard gate, the corresponding transformation is:
The transformation of all the Pauli strings are as follows:
I -----> I
Z -----> X
X -----> Z
Y -----> Y
For Phase gate, the corresponding transformation is:
The transformation of all the Pauli strings are as follows:
I -----> I
Z -----> Z
X -----> Y
Y -----> X
For CNOT gate, the corresponding transformation is:
The transformation of all the Pauli strings are as follows:
II -----> II
IZ -----> ZZ
ZI -----> ZI
ZZ -----> IZ
IX -----> IX
IY -----> ZY
ZX -----> ZX
ZY -----> IY
XI -----> XX
XZ -----> YY
YI -----> YX
YZ -----> XY
XX -----> XI
XY -----> YZ
YX -----> YI
YY -----> XZ
```

由于每个Clifford conjugation事实上都是Pauli群的一个群同构，因此箭头可以打双向，总结如下：
- Hadamard: $$ I \longleftrightarrow I \newline
Y \longleftrightarrow Y \newline
Z \longleftrightarrow X $$
- Phase: $$ I\longleftrightarrow I \newline
Z \longleftrightarrow Z \newline
X \longleftrightarrow Y $$
- CNOT: $$ II \longleftrightarrow II \newline
ZI \longleftrightarrow ZI \newline
IX \longleftrightarrow IX \newline
ZX \longleftrightarrow ZX \newline
IZ \longleftrightarrow ZZ \newline
ZY \longleftrightarrow IY \newline
XI \longleftrightarrow XX \newline
XZ \longleftrightarrow YY \newline
YI \longleftrightarrow YX \newline
YZ \longleftrightarrow XY $$

#### Result 1

考虑这样的Clifford gate：当n=2k时，取$\prod_{i=1}^k\operatorname{CNOT}_{2i-1,2i}$，当$n=2k+1$时，取$\operatorname{CNOT}_{2k,2k+1}\prod_{i=1}^k\operatorname{CNOT}_{2i-1,2i}$，发现$\mathcal{P}$的平均长度与n有简单的线性关系：
```
For 2 qubits, the average length of the transformation by [('CNOT', (0, 1))] is 1.5555555555555556
For 3 qubits, the average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (1, 2))] is 2.3333333333333335
For 4 qubits, the average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (2, 3))] is 3.111111111111111
For 5 qubits, the average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (2, 3)), ('CNOT', (3, 4))] is 3.888888888888889
For 6 qubits, the average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (2, 3)), ('CNOT', (4, 5))] is 4.666666666666667
For 7 qubits, the average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (2, 3)), ('CNOT', (4, 5)), ('CNOT', (5, 6))] is 5.444444444444445
```


#### Result 2

对于3比特情形，Result 1给出的方案并不是最好的结果，因为存在一个更好的方案如下：
```
For the Pauli operators which has initial length of 3, one can derive the following transformation:
ZZZ -----> IIZ
ZZX -----> XZX
ZZY -----> XIY
ZXZ -----> YYX
ZYZ -----> XXY
ZXX -----> ZXZ
ZXY -----> ZYI
ZYX -----> IYI
ZYY -----> IXZ
XZZ -----> ZXX
YZZ -----> IXY
XZX -----> YYZ
XZY -----> YXI
YZX -----> XYI
YZY -----> XXZ
XXZ -----> XZZ
XYZ -----> YII
YXZ -----> YZI
YYZ -----> XIZ
XXX -----> IIX
XXY -----> IZY
XYX -----> ZZY
XYY -----> ZIX
YXX -----> ZIY
YXY -----> ZZX
YYX -----> IZX
YYY -----> IIY
The average length of the transformation by [('CNOT', (0, 1)), ('CNOT', (1, 2)), ('CNOT', (2, 0))] is 2.185185185185185
```

## 代码注释（copilot生成，偷个懒……）

### 文件结构
- `Stabilizer string/pauli_string.py`

### 导入的库
```python
import numpy as np
import itertools
import matplotlib.pyplot as plt
```

### 函数和类

#### 函数
1. **`Binary_Addition_Mod2(x, y)`**
   - 功能：对两个二进制字符串进行逐位加法并取模2。
   - 示例：`Binary_Addition_Mod2("1100", "1010")` 返回 `"0110"`。

2. **`Odot(x, y, n)`**
   - 功能：检查两个字符串的反/对易关系。
   - 示例：`Odot("1100", "1010", 2)` 返回 [`1`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A50%2C%22character%22%3A6%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")。

3. **[`Generate_Pauli_Group(n)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A18%2C%22character%22%3A4%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")**
   - 功能：生成长度为2n的所有二进制字符串（n量子比特上的整个泡利群）。
   - 示例：[`Generate_Pauli_Group(2)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A18%2C%22character%22%3A4%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition") 返回 `["0000", "0001", ..., "1111"]`。

4. **[`Get_Minimun_Cases(n)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A24%2C%22character%22%3A4%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")**
   - 功能：获取n量子比特的最小平均长度（遍历整个泡利群以获取生成元），并将结果存储在列表中。

5. **`Get_Random_Generators(Z_1, pauli_group)`**
   - 功能：获取n量子比特泡利群的生成元。

6. **`H_Effect(sympletic_form, qubit)`**
   - 功能：对泡利算符应用Hadamard门。

7. **`CNOT_Effect(sympletic_form, control, target)`**
   - 功能：对泡利算符应用CNOT门。

8. **`S_Effect(sympletic_form, qubit)`**
   - 功能：对泡利算符应用S门。

9. **`Gate_Effect(gate, generators)`**
   - 功能：对克利福德群的生成元应用门操作。

10. **[`Get_CNOT_Tensor_Product(n)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A302%2C%22character%22%3A4%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")**
    - 功能：获取n量子比特上CNOT门的张量积。

#### 类
1. **[`Pauli_Operator`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A158%2C%22character%22%3A6%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")**
   - 功能：表示泡利算符。
   - 方法：
     - [`__init__(self, sympletic_form)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A159%2C%22character%22%3A8%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 初始化泡利算符。
     - [`Get_Pauli(self)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A163%2C%22character%22%3A33%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 从辛形式获取泡利算符。
     - [`Count_Length(self)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A165%2C%22character%22%3A27%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 计算泡利字符串的长度。

2. **[`Clifford_Operator`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A36%2C%22character%22%3A23%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition")**
   - 功能：表示克利福德算符。
   - 方法：
     - [`__init__(self, generators, gate_list=None, n=None)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A159%2C%22character%22%3A8%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 初始化克利福德算符。
     - [`Get_Initial_Generators(self)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A203%2C%22character%22%3A39%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 获取克利福德群的初始生成元。
     - [`Get_All_Pauli_Strings(self)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A206%2C%22character%22%3A17%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 获取克利福德算符对所有泡利字符串的作用。
     - [`Print_Clifford_Map(self, initial_length=None)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A258%2C%22character%22%3A8%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 打印克利福德算符对所有泡利字符串的作用。
     - [`Get_Average_Length(self, initial_length)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A37%2C%22character%22%3A38%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 获取克利福德算符的平均长度。
     - [`Get_Generators_From_Gate(self)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A218%2C%22character%22%3A35%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 从门操作获取克利福德群的生成元。
     - [`Plot_Distribution(self, initial_length)`](command:_github.copilot.openSymbolFromReferences?%5B%22%22%2C%5B%7B%22uri%22%3A%7B%22scheme%22%3A%22file%22%2C%22authority%22%3A%22%22%2C%22path%22%3A%22%2Fe%3A%2Fcode%20exercise%2FStabilizer%20string%2Fpauli_string.py%22%2C%22query%22%3A%22%22%2C%22fragment%22%3A%22%22%7D%2C%22pos%22%3A%7B%22line%22%3A287%2C%22character%22%3A8%7D%7D%5D%2C%2234a46f2c-d9f7-4aa7-a960-780eea62c269%22%5D "Go to definition"): 绘制泡利字符串长度的分布图。

### 示例代码
```python
if __name__ == "__main__":
    Hadamard = Clifford_Operator(None, [("H", 0)], 1)
    print("For Hadamard gate, the corresponding transformation is:")
    Hadamard.Print_Clifford_Map()

    Phase = Clifford_Operator(None, [("S", 0)], 1)
    print("For Phase gate, the corresponding transformation is:")
    Phase.Print_Clifford_Map()

    CNOT = Clifford_Operator(None, [("CNOT", (0, 1))], 2)
    print("For CNOT gate, the corresponding transformation is:")
    CNOT.Print_Clifford_Map()

    for n in range(2, 8):
        CNOT_tensor = Get_CNOT_Tensor_Product(n)
        CNOT_for_result_1 = Clifford_Operator(None, CNOT_tensor, n)
        print(f"For {n} qubits, the average length of the transformation by {CNOT_tensor} is {CNOT_for_result_1.Get_Average_Length(n)}")

    gate_list = [("CNOT", (0, 1)), ("CNOT", (1, 2)), ("CNOT", (2, 0))]
    number_of_qubits = 3
    CNOT_13_23_12 = Clifford_Operator(None, gate_list, number_of_qubits)
    CNOT_13_23_12.Print_Clifford_Map(number_of_qubits)
    print(f"The average length of the transformation by {gate_list} is {CNOT_13_23_12.Get_Average_Length(number_of_qubits)}")
    CNOT_13_23_12.Plot_Distribution(number_of_qubits)
```

### 总结
你的代码实现了对泡利算符和克利福德算符的各种操作，包括生成泡利群、应用量子门、计算平均长度和绘制分布图。这些功能对于量子计算中的稳定子码和量子纠错码的研究非常有用。
