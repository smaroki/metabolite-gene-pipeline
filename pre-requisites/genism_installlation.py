import sys
import subprocess

print("Fixing gensim installation...")
commands = [
    [sys.executable, "-m", "pip", "uninstall", "gensim", "-y"],
    [sys.executable, "-m", "pip", "install", "gensim==4.3.2"],  # Known compatible version
]

for cmd in commands:
    print(f"\nRunning: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("Error:", result.stderr)

print("\n✅ Done! Restart your kernel and try again.")
