import os, shutil, subprocess

template = open('template.vesta').read()

for dirpath, dirnames, filenames in os.walk('.'):
    for i in range(1, 10):
        stem   = f"hematite-WFN_02400_2-1_{i:01d}"
        cube   = os.path.join(dirpath, f"{stem}.cube")
        if not os.path.isfile(cube):
            print("missing", cube);  continue

        # Create a temporary .vesta scene

        scenefile = os.path.join(dirpath, f"{stem}.vesta")
        with open(scenefile, 'w') as f:
            f.write(template.replace('@@CUBE@@', cube.replace('\\','/')))

        # Render with VESTA (3.6 or newer)
        pngfile  = os.path.join(dirpath, f"{stem}.png")
        cmd = ['vesta', '-i', scenefile,
                         '-o', pngfile,
                         '-png',  # or -tiff / -jpeg

                         '-sx', '1920', '-sy', '1080']
        subprocess.run(cmd, check=True)

        os.remove(scenefile)  # tidy up