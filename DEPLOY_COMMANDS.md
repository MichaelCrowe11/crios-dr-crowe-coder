# 🚀 Deploy CriOS - Copy & Paste Commands

**Choose your deployment method and run the commands below.**

---

## ⚡ FASTEST: One-Line Deployment Script

```bash
./deploy.sh
```

This interactive script will:
- ✅ Login to Vercel (if needed)
- ✅ Deploy Frontend and/or UI
- ✅ Configure environment variables
- ✅ Optionally deploy backend to Fly.io

---

## 🎯 OPTION 1: Deploy Everything (Recommended)

### Step 1: Login to Services
```bash
# Login to Vercel
vercel login

# Login to Fly.io
curl -L https://fly.io/install.sh | sh
flyctl auth login
```

### Step 2: Deploy Frontend to Vercel
```bash
cd platform/frontend

# Set environment variables
cat > .env.production <<EOF
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev
NODE_ENV=production
EOF

# Deploy
vercel --prod
cd ../..
```

### Step 3: Deploy UI to Vercel
```bash
cd ui

# Set environment variables
cat > .env.production <<EOF
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NODE_ENV=production
EOF

# Deploy
vercel --prod
cd ..
```

### Step 4: Deploy Backend to Fly.io
```bash
cd platform/backend

# Create PostgreSQL database
flyctl postgres create --name crios-db --region iad --vm-size shared-cpu-1x --volume-size 10

# Save the DATABASE_URL from the output above, then:
# Create Redis
flyctl redis create --name crios-redis --region iad --plan free

# Save the REDIS_URL from the output above, then:
# Create app
flyctl apps create crios-backend

# Set secrets (replace with your actual values)
flyctl secrets set \
  DATABASE_URL="postgres://user:pass@host.fly.dev:5432/dbname" \
  REDIS_URL="redis://default:pass@host.fly.dev:6379" \
  ANTHROPIC_API_KEY="your-anthropic-api-key-here"

# Deploy
flyctl deploy

cd ../..
```

---

## 🌐 OPTION 2: GitHub Integration (Zero-Command Deploy)

### For Frontend (Vercel)

1. **Go to Vercel**: https://vercel.com/new
2. **Import Repository**: Click "Import Git Repository"
3. **Select**: `MichaelCrowe11/crios-dr-crowe-coder`

4. **Configure Project 1 - Frontend**:
   ```
   Project Name: crios-frontend
   Root Directory: platform/frontend
   Framework Preset: Next.js
   ```

5. **Add Environment Variables**:
   ```
   NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
   NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev
   NODE_ENV=production
   ```

6. **Click Deploy**

7. **Repeat for Project 2 - UI**:
   ```
   Project Name: crios-ui
   Root Directory: ui
   Framework Preset: Next.js

   Environment Variables:
   NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
   NODE_ENV=production
   ```

8. **Click Deploy**

### For Backend (Fly.io)

Follow Step 4 from Option 1 above, or set up GitHub Actions (see below).

---

## 🤖 OPTION 3: Automated CI/CD with GitHub Actions

### Create `.github/workflows/deploy.yml`:

```bash
mkdir -p .github/workflows

cat > .github/workflows/deploy.yml <<'EOF'
name: Deploy to Vercel and Fly.io

on:
  push:
    branches: [main]

jobs:
  deploy-backend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: superfly/flyctl-actions/setup-flyctl@master
      - run: flyctl deploy --remote-only -c platform/backend/fly.toml
        env:
          FLY_API_TOKEN: ${{ secrets.FLY_API_TOKEN }}

  deploy-frontend:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: amondnet/vercel-action@v20
        with:
          vercel-token: ${{ secrets.VERCEL_TOKEN }}
          vercel-org-id: ${{ secrets.VERCEL_ORG_ID }}
          vercel-project-id: ${{ secrets.VERCEL_FRONTEND_PROJECT_ID }}
          working-directory: ./platform/frontend
          vercel-args: '--prod'

  deploy-ui:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: amondnet/vercel-action@v20
        with:
          vercel-token: ${{ secrets.VERCEL_TOKEN }}
          vercel-org-id: ${{ secrets.VERCEL_ORG_ID }}
          vercel-project-id: ${{ secrets.VERCEL_UI_PROJECT_ID }}
          working-directory: ./ui
          vercel-args: '--prod'
EOF

# Commit and push
git add .github/workflows/deploy.yml
git commit -m "Add CI/CD workflow"
git push
```

**Required GitHub Secrets**:
- `FLY_API_TOKEN` - Get from: `flyctl auth token`
- `VERCEL_TOKEN` - Get from: https://vercel.com/account/tokens
- `VERCEL_ORG_ID` - Get from: Vercel project settings
- `VERCEL_FRONTEND_PROJECT_ID` - Get from: Vercel project settings
- `VERCEL_UI_PROJECT_ID` - Get from: Vercel project settings

---

## 🔧 Quick Commands Reference

### Vercel Commands
```bash
vercel login                    # Login to Vercel
vercel                          # Deploy to preview
vercel --prod                   # Deploy to production
vercel ls                       # List deployments
vercel logs                     # View logs
vercel env ls                   # List environment variables
vercel domains                  # Manage domains
```

### Fly.io Commands
```bash
flyctl auth login               # Login to Fly.io
flyctl deploy                   # Deploy app
flyctl status                   # Check app status
flyctl logs                     # View logs
flyctl ssh console              # SSH into app
flyctl scale count 2            # Scale to 2 instances
flyctl scale memory 2048        # Scale to 2GB RAM
```

---

## 🔑 Environment Variables

### Frontend (Vercel)
```bash
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NEXT_PUBLIC_WS_URL=wss://crios-backend.fly.dev
NODE_ENV=production
```

### UI (Vercel)
```bash
NEXT_PUBLIC_API_URL=https://crios-backend.fly.dev
NODE_ENV=production
```

### Backend (Fly.io)
```bash
DATABASE_URL=postgres://user:pass@host.fly.dev:5432/dbname
REDIS_URL=redis://default:pass@host.fly.dev:6379
ANTHROPIC_API_KEY=your-anthropic-api-key
PYTHONPATH=/app
ENV=production
```

---

## ✅ Deployment Verification

After deploying, verify everything works:

```bash
# Check frontend
curl https://crios-frontend.vercel.app

# Check UI
curl https://crios-ui.vercel.app

# Check backend health
curl https://crios-backend.fly.dev/health

# Check API docs
open https://crios-backend.fly.dev/docs
```

---

## 🆘 Troubleshooting

### "Not logged in" Error
```bash
# For Vercel
vercel login

# For Fly.io
flyctl auth login
```

### Build Fails
```bash
# Check build locally
cd platform/frontend
npm install
npm run build
```

### Backend Connection Issues
```bash
# Check backend status
flyctl status -a crios-backend

# View backend logs
flyctl logs -a crios-backend

# SSH into backend
flyctl ssh console -a crios-backend
```

### Environment Variables Not Working
```bash
# For Vercel - redeploy after setting variables
vercel --prod

# For Fly.io - verify secrets
flyctl secrets list -a crios-backend
```

---

## 📊 Expected Deployment Times

| Component | Time | Method |
|-----------|------|--------|
| Frontend (Vercel) | 2-3 min | CLI or GitHub |
| UI (Vercel) | 2-3 min | CLI or GitHub |
| Backend (Fly.io) | 5-7 min | CLI |
| Total (All) | **10-15 min** | First time |

---

## 💡 Pro Tips

1. **Use GitHub Integration**: Set it up once, auto-deploy on every push
2. **Test Locally First**: Run `npm run build` before deploying
3. **Monitor Deployments**: Check Vercel and Fly.io dashboards
4. **Set Up Alerts**: Configure monitoring for production
5. **Use Preview Deployments**: Test PRs before merging

---

## 🎯 Recommended Deployment Flow

```bash
# 1. Use the automated script
./deploy.sh

# 2. Choose "4" (Everything)
# 3. Follow prompts
# 4. Done! ✅
```

---

## 📞 Need Help?

- **Vercel Docs**: https://vercel.com/docs
- **Fly.io Docs**: https://fly.io/docs
- **Full Guide**: See `DEPLOYMENT.md`
- **Quick Guide**: See `DEPLOY_NOW.md`

---

**Ready to deploy? Run:** `./deploy.sh` 🚀
