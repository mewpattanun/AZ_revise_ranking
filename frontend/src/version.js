/* global __APP_VERSION__, __GIT_COMMIT__, __BUILD_DATE__ */
// Filled in at build time by vite.config.js (see "Build info" there).

export const APP_VERSION = __APP_VERSION__;   // from frontend/package.json, e.g. "3.0.0"
export const GIT_COMMIT = __GIT_COMMIT__;     // short git hash, "-dirty" if built with uncommitted changes
export const BUILD_DATE = __BUILD_DATE__;     // YYYY-MM-DD

// Short form for the UI, e.g. "v3.0.0"
export const VERSION_LABEL = `v${APP_VERSION}`;
