name: Build Disting NT Plugin

on:
  push:
  workflow_dispatch:

jobs:
  build:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout your repo
        uses: actions/checkout@v4

      - name: Install ARM toolchain
        run: |
          sudo apt-get update
          sudo apt-get install -y gcc-arm-none-eabi make

      - name: Clone distingNT_API
        run: |
          git clone --depth 1 https://github.com/expertsleepersltd/distingNT_API.git

      - name: Build plugin
        run: |
          make DISTINGNT_API=$PWD/distingNT_API

      - name: Upload plugin .o
        uses: actions/upload-artifact@v4
        with:
          name: distingnt-plugin
          path: |
            plugins/*.o
