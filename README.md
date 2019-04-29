# Ejemplos de Optimización Topológica en MATLAB/Octave

### Basados en el código del artículo: "Efficient topology optimization in MATLAB using 88 lines of code"

#### Elementos finitos para calor y elasticidad en 2D
Solución numérica de la ecuación de calor y elasticidad mediante elementos finitos.
* [`FEM/heat2D`](FEM/heat2D.m) 
* [`FEM/elasticity2D`](FEM/elasticity2D.m)

#### Ejemplos de Optimización
1. Ejemplo MBB  
...  
[`MBB_complete.m`](MBB_complete.m)

2. MBB simétrico  
...  
[`MBB.m`](MBB.m) (Varios ejercicios)

   * Independencia del mallado, SIMP
    
   * Patrón ajedrezado y zonas grises, Filtros

3. Diseño sencillo de una bicicleta  
...  
[`bicycle_single_load.m`](bicycle_single_load.m)

4. Ejemplo de un puente  
...  
[`bridge.m`](bridge.m)

5. Soporte en L  
...  
[`Lbracket.m`](Lbracket.m)

6. Inversor de una fuerza  
...  
[`inverter.m`](inverter.m)

7. Soporte bajo múltiples cargas  
...  
[`multiple_loads.m`](multiple_loads.m)

8. MBB robusto  
...  
[`robust_MBB.m`](robust_MBB.m)

9. Simulación para la ecuación del calor  
...  
[[images/heat_sink.png]]
[`heat_sink.m`](heat_sink.m)

10. Simulación para la ecuación del calor, con un forzante sobre un cuadrado fijo  
...  
[`square.m`](square.m)

11. Optimización en 3 dimensiones, tomado de [top3D](http://www.top3dapp.com/)  
...  
[`top3d.m`](top3d.m)

#### Archivos de MATLAB necesarios para el método MMA
Códigos en MATLAB para el método de optimización MMA (Mixed Moving Asymptotes)  
* [mmasub.m](https://pastebin.ubuntu.com/p/YNc4sg5ckB/)
* [subsolv.m](https://pastebin.ubuntu.com/p/y4pydcMWxX/)

### Two Levels Preconditioners (in progress)

1. Elasticity basis

2. Heat basis

3. [Random heat basis](Two%20Levels%20Preconditioners/Random%20heat%20basis)
