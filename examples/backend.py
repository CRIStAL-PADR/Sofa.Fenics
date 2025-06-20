import sys
from setuptools import Distribution, Extension
from setuptools.command.build_ext import build_ext
import os

def compile(name):
    # Création de l'extension C
    ext_modules = [
        Extension(
            name="fenics_tmp",
            sources=[name],
            extra_compile_args=["-O2"]
        )
    ]

    realname = name.split(".")[0]
    # Configuration de la distribution (sans setup.py)
    dist = Distribution({
        "name": "fenics_tmp",
        "ext_modules": ext_modules,
    })

    # Configuration et exécution de la commande build_ext
    cmd = build_ext(dist)
    cmd.ensure_finalized()
    cmd.run()

    # Affichage du chemin du fichier compilé
    for ext in cmd.get_outputs():
        print(f"✅ Fichier compilé : {ext}")

    import shutil
    for output in cmd.get_outputs():
        if output.endswith((".so", ".pyd", ".dll", ".dylib")):
            tokens = output.split(".")
            shortname, extension = ext_modules[0].name, tokens[-1]
            newname = "./"+shortname+"."+extension
            print(f"NAME: {newname}")
            shutil.copy(output, newname)
            print(f"📦 Copié dans ./ : {newname}")    

compile("fenics.c")