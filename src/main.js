/**
 * HTA Meta-Analysis Suite Pro
 * Main Entry Point
 */

import './styles/main.css';
import { App } from './app/app.js';

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
  App.init();
});

// Export App for console debugging
window.App = App;
