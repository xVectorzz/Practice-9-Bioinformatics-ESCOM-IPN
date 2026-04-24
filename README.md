# Práctica 9: Alineamiento Estructural y RMSD

Compara dos estructuras PDB mediante alineamiento global de secuencias para identificar pares de residuos y calcular el RMSD estructural.

## Instrucciones de uso

Asegúrate de tener `Practica9.py`, `1CHO.pdb` y `model01.pdb` en tu carpeta. Abre tu terminal y corre estos comandos en orden para preparar el entorno, instalar lo necesario y ejecutar el script:

**En Ubuntu (WSL / Linux):**
```bash
# Crear entorno, activar, instalar dependencias y ejecutar
python3 -m venv bioenv
source bioenv/bin/activate
pip install biopython
python3 Practica9.py
```

**En Windows (CMD):**
```cmd
# Crear entorno, activar, instalar dependencias y ejecutar
python -m venv bioenv
bioenv\Scripts\activate
pip install biopython
python Practica9.py
```
