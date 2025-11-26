import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['src/**/__tests__/*.test.js'],
    coverage: {
      reporter: ['text', 'html'],
      include: ['src/engine/**/*.js'],
      exclude: ['src/engine/__tests__/**']
    }
  }
});
