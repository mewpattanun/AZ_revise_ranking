import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import { execSync } from "child_process";
import { readFileSync } from "fs";
import { resolve } from "path";

// ── Build info shown in the header badge ─────────────────────
// Version: bump "version" in package.json when you release a change
// (3.0.1 = fix, 3.1.0 = new feature or scoring change, 4.0.0 = big redesign).
// Commit and date are filled in automatically every time the site is built,
// so the live site always says exactly which code it is running.
const pkg = JSON.parse(readFileSync(resolve(__dirname, "package.json"), "utf-8"));

function git(args) {
  try {
    return execSync(`git ${args}`, { cwd: __dirname, stdio: ["ignore", "pipe", "ignore"] })
      .toString()
      .trim();
  } catch {
    return ""; // not a git checkout (e.g. a zip download)
  }
}

// "-dirty" = built from code with uncommitted changes. Files that the deploy
// itself rewrites (package-lock.json after npm install, the dist build output)
// are ignored, so a clean deploy on the server never shows "-dirty".
const changed = git("status --porcelain --untracked-files=no")
  .split("\n")
  .filter((line) => line && !/frontend\/(package-lock\.json|dist\/)/.test(line));
const commit = (git("rev-parse --short HEAD") || "unknown") + (changed.length ? "-dirty" : "");

export default defineConfig({
  plugins: [react()],
  define: {
    __APP_VERSION__: JSON.stringify(pkg.version),
    __GIT_COMMIT__: JSON.stringify(commit),
    __BUILD_DATE__: JSON.stringify(new Date().toISOString().slice(0, 10)),
  },
  build: {
    rollupOptions: {
      input: {
        main: resolve(__dirname, "index.html"),
        test: resolve(__dirname, "test.html"),
      },
    },
  },
  server: {
    port: 3000,
    proxy: {
      "/api": {
        target: "http://localhost:8000",
        changeOrigin: true,
        timeout: 300000,
        configure: (proxy) => {
          proxy.on("proxyReq", (proxyReq) => {
            proxyReq.socket?.setTimeout(300000);
          });
        },
      },
    },
  },
});
