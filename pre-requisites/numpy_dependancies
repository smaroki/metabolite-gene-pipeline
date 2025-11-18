import sys
import subprocess

print("🔧 Fixing numpy and dependencies...")

# Uninstall conflicting packages
packages_to_remove = ["numpy", "scipy", "gensim", "pandas"]
for pkg in packages_to_remove:
    print(f"\nUninstalling {pkg}...")
    subprocess.run([sys.executable, "-m", "pip", "uninstall", pkg, "-y"], 
                   capture_output=True)

# Reinstall with compatible versions
print("\n📦 Installing compatible versions...")
packages_to_install = [
    "numpy==1.24.3",
    "scipy==1.11.4", 
    "pandas==2.0.3",
    "gensim==4.3.2"
]

for pkg in packages_to_install:
    print(f"\nInstalling {pkg}...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", pkg, "--force-reinstall", "--no-cache-dir"],
        capture_output=True, 
        text=True
    )
    if "Successfully installed" in result.stdout:
        print(f"  ✅ {pkg} installed")
    else:
        print(f"  ⚠️ {result.stdout}")

print("\n✅ Done! Now RESTART YOUR KERNEL!")
print("   Kernel → Restart Kernel")
