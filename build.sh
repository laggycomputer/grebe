cargo build --release

other_targets=(
    # legacy raspberry pi, the most comical of benchmarking systems
    # sudo apt install gcc-arm-linux-gnueabihf
    "armv7-unknown-linux-gnueabihf"
    # sudo apt install mingw-w64
    "x86_64-pc-windows-gnu"
    "i686-pc-windows-gnu"
    # i love apple's ARR policies
#     "x86_64-apple-darwin"
)
for target in "${other_targets[@]}"; do
    rustup target add $target
    cargo build --target=$target --release
done
