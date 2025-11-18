import sys
import subprocess

print("🧹 STEP 1: Complete cleanup...")
print("=" * 60)

# Remove ALL conflicting packages
packages_to_remove = [
    "pandas", "numpy", "scipy", "gensim", "ms2query", 
    "matchms", "scikit-learn", "matplotlib"
]

for pkg in packages_to_remove:
    print(f"Removing {pkg}...", end=" ")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "uninstall", pkg, "-y"],
        capture_output=True
    )
    print("✓")

print("\n📦 STEP 2: Installing compatible versions...")
print("=" * 60)

# Install in specific order with compatible versions
install_order = [
    "numpy==1.24.3",
    "scipy==1.10.1", 
    "pandas==2.0.3",
    "gensim==4.3.2",
    "scikit-learn==1.3.2",
    "matchms==0.26.4",
    "ms2query==1.5.4"
]

for pkg in install_order:
    print(f"\nInstalling {pkg}...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", pkg, "--no-cache-dir"],
        capture_output=True,
        text=True
    )
    
    if "Successfully installed" in result.stdout or "Requirement already satisfied" in result.stdout:
        print(f"  ✅ {pkg.split('==')[0]} installed")
    else:
        print(f"  ❌ Error: {result.stderr[:100]}")

print("\n" + "=" * 60)
print("✅ INSTALLATION COMPLETE!")
print("=" * 60)
print("\n⚠️  CRITICAL: You MUST restart your kernel now!")
print("   Click: Kernel → Restart Kernel")
