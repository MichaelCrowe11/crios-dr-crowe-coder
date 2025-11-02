#!/bin/bash
# CriOS Vercel Deployment Script
# This script helps deploy the frontend components to Vercel

set -e

echo "🚀 CriOS Vercel Deployment Script"
echo "=================================="
echo ""

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo -e "${YELLOW}⚠️  Vercel CLI not found. Installing...${NC}"
    npm i -g vercel
fi

# Menu
echo "Select deployment target:"
echo "1) Frontend (platform/frontend)"
echo "2) Algorithm Studio UI (ui)"
echo "3) Both"
echo ""
read -p "Enter choice [1-3]: " choice

deploy_frontend() {
    echo ""
    echo -e "${BLUE}📦 Deploying Frontend...${NC}"
    cd platform/frontend

    read -p "Deploy to production? (y/n): " prod
    if [ "$prod" = "y" ] || [ "$prod" = "Y" ]; then
        vercel --prod
    else
        vercel
    fi

    cd ../..
    echo -e "${GREEN}✅ Frontend deployed!${NC}"
}

deploy_ui() {
    echo ""
    echo -e "${BLUE}📦 Deploying Algorithm Studio UI...${NC}"
    cd ui

    read -p "Deploy to production? (y/n): " prod
    if [ "$prod" = "y" ] || [ "$prod" = "Y" ]; then
        vercel --prod
    else
        vercel
    fi

    cd ..
    echo -e "${GREEN}✅ UI deployed!${NC}"
}

case $choice in
    1)
        deploy_frontend
        ;;
    2)
        deploy_ui
        ;;
    3)
        deploy_frontend
        deploy_ui
        ;;
    *)
        echo -e "${YELLOW}Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

echo ""
echo -e "${GREEN}🎉 Deployment complete!${NC}"
echo ""
echo "Next steps:"
echo "1. Set environment variables in Vercel Dashboard"
echo "2. Deploy backend to Fly.io (see DEPLOYMENT.md)"
echo "3. Configure CORS on backend"
echo "4. Test the deployment"
echo ""
echo "URLs:"
echo "  Frontend: https://crios-frontend.vercel.app"
echo "  UI: https://crios-ui.vercel.app"
echo "  Dashboard: https://vercel.com/dashboard"
