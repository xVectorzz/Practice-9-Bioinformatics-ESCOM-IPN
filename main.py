#!/usr/bin/env python3
"""
auto_select_and_rmsd.py

Detecta la mejor correspondencia entre cadenas de un PDB de referencia
y un modelo (dock), muestra coincidencias por % identidad, y calcula RMSD
sobre átomos C-alpha emparejados por alineamiento.
"""

from Bio.PDB import PDBParser, Superimposer, is_aa
from Bio import pairwise2
from Bio.Data.IUPACData import protein_letters_3to1
import sys

# ---------- CONFIG ----------
pdb_ref = "1CHO.pdb"            # tu PDB de referencia (ej. 1CHO)
pdb_model = "model01.pdb"  # tu modelo/docking (ej. ClusPro)
# Si quieres forzar cadenas, pon aquí (ej "E","A"). Dejar None para auto
force_chain_ref = None
force_chain_model = None
# -----------------------------

parser = PDBParser(QUIET=True)

# --- utilidades ---
def three_to_one(resname):
    """Convierte nombre de 3 letras a 1 letra; devuelve 'X' si desconocido."""
    try:
        return protein_letters_3to1[resname.capitalize()]
    except KeyError:
        return "X"

def ca_residue_list(chain):
    """Devuelve lista de (residue_obj, 1-letter AA) solo para residuos con CA y que sean AA."""
    out = []
    for res in chain:
        # res.id[0] == ' ' means standard residue (no hetero)
        if not is_aa(res, standard=True):
            continue
        if 'CA' in res:
            one = three_to_one(res.get_resname())
            out.append((res, one))
    return out

def seq_from_ca_list(ca_list):
    return "".join([x[1] for x in ca_list])

# --- cargar estructuras ---
try:
    struct_ref = parser.get_structure("ref", pdb_ref)[0]
    struct_model = parser.get_structure("model", pdb_model)[0]
except FileNotFoundError as e:
    print("Error: no se encontró un archivo PDB. Revisa las rutas:", e)
    sys.exit(1)

# listar cadenas
chains_ref = list(struct_ref.get_chains())
chains_model = list(struct_model.get_chains())

print("Cadenas en referencia:", [c.id for c in chains_ref])
print("Cadenas en modelo  :", [c.id for c in chains_model])
print()

# si forzamos cadenas, valida y selecciona
if force_chain_ref:
    try:
        chains_ref = [struct_ref[force_chain_ref]]
    except KeyError:
        print(f"Cadena {force_chain_ref} no encontrada en {pdb_ref}")
        sys.exit(1)
if force_chain_model:
    try:
        chains_model = [struct_model[force_chain_model]]
    except KeyError:
        print(f"Cadena {force_chain_model} no encontrada en {pdb_model}")
        sys.exit(1)

# Construir secuencias CA por cadena
chain_info_ref = []
for c in chains_ref:
    ca_list = ca_residue_list(c)
    seq = seq_from_ca_list(ca_list)
    chain_info_ref.append((c.id, c, ca_list, seq))

chain_info_model = []
for c in chains_model:
    ca_list = ca_residue_list(c)
    seq = seq_from_ca_list(ca_list)
    chain_info_model.append((c.id, c, ca_list, seq))

# Comparar todas las parejas y calcular % identidad basado en alineamiento global
results = []
for id_r, chain_r, ca_r, seq_r in chain_info_ref:
    for id_m, chain_m, ca_m, seq_m in chain_info_model:
        if len(seq_r) == 0 or len(seq_m) == 0:
            pid = 0.0
            aln_len = 0
            score = 0
        else:
            # alineamiento global simple (xx cuenta matches)
            aln = pairwise2.align.globalxx(seq_r, seq_m, one_alignment_only=True)[0]
            a_r, a_m = aln.seqA, aln.seqB
            matches = sum(1 for a,b in zip(a_r, a_m) if a == b and a != '-' and b != '-')
            aligned_positions = sum(1 for a,b in zip(a_r, a_m) if a != '-' and b != '-')
            pid = (matches / aligned_positions)*100 if aligned_positions>0 else 0.0
            aln_len = aligned_positions
            score = aln.score
        results.append({
            "chain_ref": id_r,
            "chain_model": id_m,
            "len_ref": len(seq_r),
            "len_model": len(seq_m),
            "aligned_pos": aln_len,
            "identity_pct": pid,
            "score": score
        })

# Ordenar por identidad descendente
results_sorted = sorted(results, key=lambda x: (x["identity_pct"], x["aligned_pos"]), reverse=True)

# Mostrar tabla de top coincidencias
print("Top coincidencias (por % identidad):")
print(f"{'Ref':>4} {'Model':>6} {'len_ref':>7} {'len_mod':>7} {'aligned':>7} {'%id':>7}")
for r in results_sorted[:10]:
    print(f"{r['chain_ref']:>4} {r['chain_model']:>6} {r['len_ref']:7d} {r['len_model']:7d} {r['aligned_pos']:7d} {r['identity_pct']:7.2f}")

# Seleccionar pareja a usar: la de mayor identidad (primera de results_sorted)
best = results_sorted[0]
sel_ref = best['chain_ref']
sel_model = best['chain_model']
print()
print("Pareja seleccionada automáticamente:", sel_ref, " (ref)  <-> ", sel_model, " (model)")
print(f"Identidad {best['identity_pct']:.2f}%  (pos alineadas: {best['aligned_pos']})")
print()

# Alternativa: permitir que el usuario elija otra pareja
choice = input("¿Usar esta pareja? (Y/n): ").strip().lower()
if choice == 'n':
    print("Introduce la pareja manual en formato REF,MODEL (ej: E,A):")
    manual = input("Pareja: ").strip()
    try:
        r_id, m_id = [x.strip() for x in manual.split(",")]
        sel_ref = r_id
        sel_model = m_id
    except Exception as e:
        print("Entrada inválida. Saliendo.")
        sys.exit(1)

# Obtener objetos chain y sus listas CA
chain_r_obj = struct_ref[sel_ref]
chain_m_obj = struct_model[sel_model]
ca_r = ca_residue_list(chain_r_obj)
ca_m = ca_residue_list(chain_m_obj)
seq_r = seq_from_ca_list(ca_r)
seq_m = seq_from_ca_list(ca_m)

# Re-obtener alineamiento para construir pares Cα
aln = pairwise2.align.globalxx(seq_r, seq_m, one_alignment_only=True)[0]
a_r, a_m = aln.seqA, aln.seqB

paired_ca_ref = []
paired_ca_model = []
i = j = 0
for ar, am in zip(a_r, a_m):
    if ar != '-' and am != '-':
        # ambos corresponden: añadir CA
        paired_ca_ref.append(ca_r[i][0]['CA'])
        paired_ca_model.append(ca_m[j][0]['CA'])
    if ar != '-':
        i += 1
    if am != '-':
        j += 1

print(f"Residues emparejados (CA): {len(paired_ca_ref)}")
if len(paired_ca_ref) < 4:
    print("Advertencia: muy pocos pares emparejados para un RMSD confiable (<4).")
    # continuar igualmente

# Calcular RMSD
sup = Superimposer()
sup.set_atoms(paired_ca_ref, paired_ca_model)
sup.apply(struct_model.get_atoms())  # aplica transform al modelo
print(f"RMSD (Cα, basado en alineamiento): {sup.rms:.3f} Å")

# Fin
