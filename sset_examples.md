

```python
import sset as ss
```


```python
Z = ss.scomplex(["xy","yz","xz","zw","xyz","yzw","xzw","xyw"])
#Z is topologically a sphere

print("Sigma_0 = ", Z.Sigma(0),"\n")
print("Sigma_1 = ", Z.Sigma(1),"\n")
print("Sigma_2 = ", Z.Sigma(2),"\n")

print("Betti numbers =",Z.betti())
```

    Sigma_0 =  ['y', 'z', 'w', 'x'] 
    
    Sigma_1 =  ['xy', 'zw', 'yw', 'yz', 'xz', 'xw'] 
    
    Sigma_2 =  ['xzw', 'xyw', 'yzw', 'xyz'] 
    
    Betti numbers = [1, 0, 1]



```python
X = ss.semi_sset([{"x","y","z"},{"xy","yz","xz"}], [["xy",["x",[]],["y",[]]],["yz",["y",[]],["z",[]]],["xz",["x",[]],["z",[]]]])
#X is a circle
print(X.betti())

T = ss.semi_sset([{'x'},{'a','b','c'},{'A','B'}],[['A',['a',['x',[]],['x',[]]],['b',['x',[]],['x',[]]],['c',['x',[]],['x',[]]]],['B', ['a',['x',[]],['x',[]]],['b',['x',[]],['x',[]]],['c',['x',[]],['x',[]]]]])
#T is a torus
print(T.betti())

```

    [1, 1]
    [1, 2, 1]



```python
X = ss.sset([{"x","y","z"},{"a","b","c"},{"A"}], ["A",["a",["x",[]],["y",[]]],["b",["y",[]],["z",[]]],["c",["x",[]],["z",[]]]])
# X is a triangle (homemorphic to a closed disc, hence contractible)
print(X.simplices)
print(X.structuretree)

print("Betti numbers =", X.betti())

Q = X.simplexcollapse("a")
print(Q.simplices)
print(Q.structuretree)

print("Betti numbers=", Q.betti())

Q1 = Q.simplexcollapse("b")

print(Q1.simplices)
print(Q1.structuretree)
print("Betti numbers=", Q1.betti())

Q2 = Q1.simplexcollapse("c")
#Q_2 is homeomorphic to a sphere

print(Q2.simplices)
print(Q2.structuretree)
print("Betti numbers =", Q2.betti())

```

    [{'y', 'x', 'z'}, {'a', 'b', 'c'}, {'A'}]
    ['A', ['a', ['x', []], ['y', []]], ['b', ['y', []], ['z', []]], ['c', ['x', []], ['z', []]]]
    Betti numbers = [1, 0, 0]
    
    collapsing simplex a....
    
    [{'a_0', 'z'}, {'b', 'c', '0-a_0'}, {'A'}]
    ['A', ['0-a_0', ['a_0', []], ['a_0', []]], ['b', ['a_0', []], ['z', []]], ['c', ['a_0', []], ['z', []]]]
    Betti numbers= [1, 0, 0]
    
    collapsing simplex b....
    
    [{'b_0'}, {'0-b_0', 'c'}, {'A'}]
    ['A', ['0-b_0', ['b_0', []], ['b_0', []]], ['0-b_0', ['b_0', []], ['b_0', []]], ['c', ['b_0', []], ['b_0', []]]]
    Betti numbers= [1, 0, 0]
    
    collapsing simplex c....
    
    [{'c_0'}, {'0-c_0'}, {'A'}]
    ['A', ['0-c_0', ['c_0', []], ['c_0', []]], ['0-c_0', ['c_0', []], ['c_0', []]], ['0-c_0', ['c_0', []], ['c_0', []]]]
    Betti numbers = [1, 0, 1]



```python

```
